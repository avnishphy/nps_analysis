#!/usr/bin/env bash
# Submit one compare_livetimes SWIF2 job per run and one dependent finalizer
# per kinematic setting. The finalizer validates/merges CSV parts, makes the
# three-page PDF, and writes an auditable result tar.
set -euo pipefail

SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
DEFAULT_ROOT="/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_final"
CONFIG_REL="config/nps_dvcs_all_kins_main.csv"
EFF_REL="src/efficiencies"
EXE_REL="${EFF_REL}/compare_livetimes"
SOURCE_REL="${EFF_REL}/compare_livetimes.cxx"
PLOT_REL="${EFF_REL}/plot_compare_livetimes.py"

die() { echo "[ERROR] $*" >&2; exit 1; }
need_value() { [[ $# -ge 2 && -n "$2" && "$2" != --* ]] || die "$1 requires a value"; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Command not found: $1"; }
safe_name() { printf '%s' "$1" | sed 's/[^[:alnum:]_-]/_/g'; }

farm_reexec() {
  [[ "${NPS_FARM_ENV_LOADED:-0}" == 1 ]] && return 0
  local reentry
  reentry="$(mktemp "${TMPDIR:-/tmp}/compare_livetimes_reentry.XXXXXX.sh")"
  {
    printf '#!/usr/bin/env bash\nset -euo pipefail\nself=%q\nrm -f "$self"\nexport NPS_FARM_ENV_LOADED=1\nexec %q' \
      "$reentry" "$SCRIPT_PATH"
    printf ' %q' "$@"
    printf '\n'
  } > "$reentry"
  chmod 700 "$reentry"
  export NPS_FARM_REENTRY="$reentry"
  exec csh -f -c 'if ( -f /usr/share/Modules/init/csh ) source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; if ( $status != 0 ) exit 97; exec /bin/bash "$NPS_FARM_REENTRY"'
}

prepare_runtime() {
  local runtime_tar="$1" root_name="$2"
  [[ "$runtime_tar" != */* ]] || die "Staged runtime tar must be a basename"
  [[ "$root_name" =~ ^[[:alnum:]_.-]+$ ]] || die "Unsafe runtime root name"
  [[ -s "$runtime_tar" && -s "${runtime_tar}.sha256" ]] || die "Runtime tar/checksum missing"
  sha256sum -c "${runtime_tar}.sha256"
  tar -tzf "$runtime_tar" >/dev/null
  tar -tzf "$runtime_tar" | grep -Fx "${root_name}/${EXE_REL}" >/dev/null || die "Runtime lacks executable"
  tar -tzf "$runtime_tar" | grep -Fx "${root_name}/${PLOT_REL}" >/dev/null || die "Runtime lacks plotter"
  tar -xzf "$runtime_tar"
}

# SWIF worker mode: staged manifests name the exact, ordered replay segments.
if [[ "${1:-}" == --run-one ]]; then
  worker_args=("$@")
  [[ $# -eq 5 ]] || die "--run-one requires TAR ROOT_NAME KIN RUN"
  runtime_tar="$2"; root_name="$3"; kin="$4"; run="$5"
  [[ "$run" =~ ^[0-9]+$ ]] || die "Invalid run: $run"
  [[ -s input_manifest.tsv && -s run_source_report.csv ]] || die "Input metadata missing"

  farm_reexec "${worker_args[@]}"
  prepare_runtime "$runtime_tar" "$root_name"
  runtime_root="${PWD}/${root_name}"
  exe="${runtime_root}/${EXE_REL}"
  config="${PWD}/compare_config.csv"
  [[ -x "$exe" && -s "$config" ]] || die "Runtime executable or staged config unavailable"

  files=(); source_kind=""; input_access=""
  while IFS=$'\t' read -r row_run row_source row_access basename source_path size_bytes; do
    [[ "$row_run" == run ]] && continue
    [[ "$row_run" == "$run" ]] || die "Manifest run mismatch: $row_run"
    [[ "$row_source" == updated || "$row_source" == production ]] || die "Invalid source: $row_source"
    [[ -z "$source_kind" || "$source_kind" == "$row_source" ]] || die "Mixed replay sources"
    [[ -z "$input_access" || "$input_access" == "$row_access" ]] || die "Mixed input access"
    source_kind="$row_source"; input_access="$row_access"
    if [[ "$row_access" == stage ]]; then file="${PWD}/${basename}"
    elif [[ "$row_access" == shared ]]; then file="$source_path"
    else die "Invalid input access: $row_access"
    fi
    [[ -r "$file" ]] || die "ROOT input unreadable: $file"
    [[ "$(stat -c %s "$file")" == "$size_bytes" ]] || die "ROOT input size changed: $file"
    files+=("$file")
  done < input_manifest.tsv
  (( ${#files[@]} > 0 )) || die "Manifest contains no ROOT files"

  mkdir -p job_output
  start_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  command=("$exe" "$run" "${files[@]}" --config "$config" --output-dir "${PWD}/job_output" --file-source "$source_kind")
  printf '[worker] command:'; printf ' %q' "${command[@]}"; printf '\n'
  "${command[@]}"
  produced="${PWD}/job_output/compare_livetimes_${kin}.csv"
  [[ -s "$produced" && "$(wc -l < "$produced")" -eq 2 ]] || die "Expected one two-line CSV: $produced"
  cp "$produced" compare_livetimes_part.csv
  {
    printf 'kin=%s\nrun=%s\nsource=%s\ninput_access=%s\ninput_files=%s\n' \
      "$kin" "$run" "$source_kind" "$input_access" "${#files[@]}"
    printf 'start_utc=%s\nend_utc=%s\nruntime_sha256=%s\ncommand=' \
      "$start_utc" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$(sha256sum "$runtime_tar" | awk '{print $1}')"
    printf '%q ' "${command[@]}"; printf '\n'
  } > compare_livetimes_job_manifest.txt
  exit 0
fi

# Finalizer mode: all antecedent run outputs already reside on /volatile.
if [[ "${1:-}" == --finalize-one ]]; then
  worker_args=("$@")
  [[ $# -eq 7 ]] || die "--finalize-one requires TAR ROOT_NAME KIN OUTPUT_BASE RESULT_TAR PARTS_DIR"
  runtime_tar="$2"; root_name="$3"; kin="$4"; output_base="$5"; result_tar="$6"; parts_dir="$7"
  [[ "$output_base" == /volatile/* && "$result_tar" == /volatile/* && "$parts_dir" == /volatile/* ]] || \
    die "Final outputs must be under /volatile"
  [[ -s run_source_report.csv && -s compare_config.csv ]] || die "Staged report/config missing"

  farm_reexec "${worker_args[@]}"
  prepare_runtime "$runtime_tar" "$root_name"
  runtime_root="${PWD}/${root_name}"
  plotter="${runtime_root}/${PLOT_REL}"
  support="${runtime_root}/output/efficiency_stuff"
  safe="$(safe_name "$kin")"
  efficiency_out="${output_base}/efficiency_stuff"
  plot_out="${output_base}/plots"
  metadata_out="${output_base}/metadata/${safe}"
  mkdir -p "$efficiency_out" "$plot_out" "$metadata_out"

  merged="${efficiency_out}/compare_livetimes_${kin}.csv"
  merged_tmp="${merged}.tmp.$$"
  python3 - "$parts_dir" run_source_report.csv "$kin" "$merged_tmp" <<'PY'
import csv, pathlib, re, sys
parts_dir, report_path, kin, output = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), sys.argv[3], pathlib.Path(sys.argv[4])
with report_path.open(newline='') as f:
    expected = {int(r['run']) for r in csv.DictReader(f) if r['kin'] == kin and r['input_status'] == 'ready'}
if not expected:
    raise SystemExit(f'No expected runs for {kin}')
rows, header = {}, None
for path in sorted(parts_dir.glob('run*.csv')):
    match = re.fullmatch(r'run([0-9]+)\.csv', path.name)
    if not match:
        continue
    run = int(match.group(1))
    with path.open(newline='') as f:
        data = list(csv.reader(f))
    if len(data) != 2:
        raise SystemExit(f'{path}: expected header plus one row')
    if header is None:
        header = data[0]
    elif data[0] != header:
        raise SystemExit(f'{path}: CSV header mismatch')
    if int(data[1][0]) != run:
        raise SystemExit(f'{path}: row run does not match filename')
    if run in rows:
        raise SystemExit(f'duplicate run {run}')
    rows[run] = data[1]
missing, extra = sorted(expected - rows.keys()), sorted(rows.keys() - expected)
if missing or extra:
    raise SystemExit(f'part mismatch: missing={missing} extra={extra}')
with output.open('w', newline='') as f:
    writer = csv.writer(f, lineterminator='\n')
    writer.writerow(header)
    writer.writerows(rows[run] for run in sorted(rows))
print(f'[merge] {output}: {len(rows)} runs')
PY
  mv "$merged_tmp" "$merged"

  for prefix in efficiency selection_report; do
    src="${support}/${prefix}_${kin}.csv"
    [[ -s "$src" ]] || die "Runtime support CSV missing: $src"
    cp "$src" "$efficiency_out/"
  done
  export MPLCONFIGDIR="${PWD}/matplotlib"
  python3 "$plotter" --output-dir "$plot_out" "$merged"
  shopt -s nullglob
  pdfs=("${plot_out}/compare_livetimes_multipanel_${safe}_"*.pdf)
  (( ${#pdfs[@]} > 0 )) || die "Plotter produced no PDF for $kin"

  cp run_source_report.csv "$metadata_out/run_source_report.csv"
  cp compare_config.csv "$metadata_out/compare_config.csv"
  {
    printf 'kin=%s\nruns=%s\npdfs=%s\nruntime_sha256=%s\ncompleted_utc=%s\n' \
      "$kin" "$(( $(wc -l < "$merged") - 1 ))" "${#pdfs[@]}" \
      "$(sha256sum "$runtime_tar" | awk '{print $1}')" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  } > "$metadata_out/finalize_manifest.txt"

  bundle="${PWD}/compare_livetimes_${safe}"
  mkdir -p "$bundle/efficiency_stuff" "$bundle/plots" "$bundle/metadata"
  cp "$merged" "${efficiency_out}/efficiency_${kin}.csv" \
     "${efficiency_out}/selection_report_${kin}.csv" "$bundle/efficiency_stuff/"
  cp "${pdfs[@]}" "$bundle/plots/"
  cp "$metadata_out/run_source_report.csv" "$metadata_out/compare_config.csv" \
     "$metadata_out/finalize_manifest.txt" "$bundle/metadata/"
  [[ -d "${output_base}/metadata/${safe}/runs" ]] && cp -a "${output_base}/metadata/${safe}/runs" "$bundle/metadata/"
  result_tmp="${result_tar}.tmp.$$"
  tar -C "$PWD" -czf "$result_tmp" "$(basename "$bundle")"
  tar -tzf "$result_tmp" >/dev/null
  mv "$result_tmp" "$result_tar"
  echo "[finalizer] merged CSV: $merged"
  echo "[finalizer] result tar: $result_tar"
  exit 0
fi

usage() {
  cat <<EOF
Usage: $(basename "$0") (--kin KIN... | --all-kins | --run RUN...) [options]

Selection/input:
  --kin NAME             Select Kin_old; repeatable, comma lists accepted
  --all-kins             Select all settings matching other filters
  --run N                Restrict runs; repeatable, comma lists accepted
  --target NAME          Restrict target; repeatable (default: all)
  --types CSV            Config Type filter (default: production,Production)
  --config FILE          Run metadata CSV
  --updated-dir DIR      Updated replay directory
  --production-dir DIR   Production replay directory
  --input-access MODE    shared (default) or stage
  --missing-input MODE   skip (default) or fail

SWIF2/resources:
  --workflow NAME        Default: compare_livetimes_<timestamp>
  --account NAME         Default: hallc
  --partition NAME       Default: production
  --cores N              Cores per run job (default: 1)
  --ram SIZE             RAM per job (default: 4g)
  --disk SIZE            Disk per job (default: 10g; increase for staged ROOT files)
  --time DURATION        Default: 12h
  --max-concurrent N     Maximum simultaneous run jobs per kin (default: 8)
  --output-dir DIR       Persistent result directory under /volatile
  --tar-dir DIR          Runtime/input-manifest directory
  --analysis-root DIR    Source tree to package
  --create-only          Create workflows without starting them
  --dry-run              Build/verify inputs and tar; print SWIF2 commands
  -h, --help             Show this help

Each kinematic gets its own workflow. Run jobs use exact preflighted replay
segments and export one CSV part. A dependent finalizer requires every part,
merges them in run order, copies efficiency/selection support CSVs, generates
the run/current/S1X-rate multipage PDF, and writes a result tar. Logs go to
/farm_out; persistent parts, merged outputs, plots, metadata and tar files go
under --output-dir.
EOF
}

ANALYSIS_ROOT="$DEFAULT_ROOT"
CONFIG=""
UPDATED_DIR="/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/updated"
PRODUCTION_DIR="/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/production"
TYPES="production,Production"
INPUT_ACCESS=shared; MISSING_INPUT=skip
WORKFLOW="compare_livetimes_$(date +%Y%m%d_%H%M%S)"
ACCOUNT=hallc; PARTITION=production
CORES=1; RAM=4g; DISK=10g; WALLTIME=12h; MAX_CONCURRENT=8
OUTPUT_DIR=""; TAR_DIR=""
ALL_KINS=0; CREATE_ONLY=0; DRY_RUN=0
KINS=(); RUNS=(); TARGETS=()

while (($#)); do
  case "$1" in
    --kin|--run|--target)
      option="$1"; need_value "$@"; IFS=',' read -r -a values <<< "$2"
      case "$option" in --kin) KINS+=("${values[@]}");; --run) RUNS+=("${values[@]}");; --target) TARGETS+=("${values[@]}");; esac
      shift 2 ;;
    --all-kins) ALL_KINS=1; shift ;;
    --types) need_value "$@"; TYPES="$2"; shift 2 ;;
    --config) need_value "$@"; CONFIG="$2"; shift 2 ;;
    --updated-dir) need_value "$@"; UPDATED_DIR="$2"; shift 2 ;;
    --production-dir) need_value "$@"; PRODUCTION_DIR="$2"; shift 2 ;;
    --input-access) need_value "$@"; INPUT_ACCESS="$2"; shift 2 ;;
    --missing-input) need_value "$@"; MISSING_INPUT="$2"; shift 2 ;;
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

(( ALL_KINS || ${#KINS[@]} || ${#RUNS[@]} )) || die "Select --kin, --all-kins, or --run"
(( ALL_KINS == 0 || ${#KINS[@]} == 0 )) || die "Use --kin or --all-kins"
[[ "$INPUT_ACCESS" == shared || "$INPUT_ACCESS" == stage ]] || die "Invalid --input-access"
[[ "$MISSING_INPUT" == skip || "$MISSING_INPUT" == fail ]] || die "Invalid --missing-input"
[[ "$WORKFLOW" =~ ^[[:alnum:]][[:alnum:]_-]*$ ]] || die "Unsafe workflow name"
for value in "$CORES" "$MAX_CONCURRENT"; do [[ "$value" =~ ^[1-9][0-9]*$ ]] || die "Cores/concurrency must be positive"; done
[[ "$RAM" =~ ^[1-9][0-9]*[kmgt]?$ && "$DISK" =~ ^[1-9][0-9]*[kmgt]?$ ]] || die "Invalid RAM/disk size"
[[ "$WALLTIME" =~ ^[1-9][0-9]*[smhd]$ ]] || die "Invalid walltime"
for run in "${RUNS[@]}"; do [[ "$run" =~ ^[0-9]+$ ]] || die "Invalid run: $run"; done

ANALYSIS_ROOT="$(readlink -f "$ANALYSIS_ROOT")"
CONFIG="$(readlink -f "${CONFIG:-${ANALYSIS_ROOT}/${CONFIG_REL}}")"
UPDATED_DIR="$(readlink -f "$UPDATED_DIR")"; PRODUCTION_DIR="$(readlink -f "$PRODUCTION_DIR")"
[[ -s "$CONFIG" ]] || die "Config missing: $CONFIG"
for file in "${ANALYSIS_ROOT}/${SOURCE_REL}" "${ANALYSIS_ROOT}/${PLOT_REL}"; do [[ -s "$file" ]] || die "Missing: $file"; done

SUBMIT_USER="${NPS_SWIF_USER:-${LOGNAME:-${USER:-}}}"
[[ "$SUBMIT_USER" =~ ^[[:alpha:]_][[:alnum:]_-]*$ ]] || die "Set NPS_SWIF_USER"
OUTPUT_DIR="${OUTPUT_DIR:-/volatile/hallc/nps/${SUBMIT_USER}/compare_livetimes/${WORKFLOW}}"
TAR_DIR="$(readlink -m "${TAR_DIR:-/group/nps/${SUBMIT_USER}/swif_inputs}")"
LOG_DIR="/farm_out/${SUBMIT_USER}/compare_livetimes/${WORKFLOW}"
[[ "$OUTPUT_DIR" == /volatile/* && "$LOG_DIR" == /farm_out/* ]] || die "Unsafe output/log path"

need_cmd python3; need_cmd tar; need_cmd sha256sum; need_cmd g++; need_cmd csh
work="$(mktemp -d)"
cleanup_files=()
cleanup() {
  rm -rf -- "$work"
  for file in "${cleanup_files[@]}"; do [[ ! -e "$file" ]] || rm -f -- "$file"; done
}
trap cleanup EXIT

kin_csv="$(IFS=,; echo "${KINS[*]}")"; run_csv="$(IFS=,; echo "${RUNS[*]}")"; target_csv="$(IFS=,; echo "${TARGETS[*]}")"
python3 - "$CONFIG" "$kin_csv" "$run_csv" "$target_csv" "$TYPES" > "$work/jobs.tsv" <<'PY'
import csv, sys
path, kin_arg, run_arg, target_arg, type_arg = sys.argv[1:]
kins, runs = set(filter(None, kin_arg.split(','))), set(filter(None, run_arg.split(',')))
targets = {x.casefold() for x in filter(None, target_arg.split(',')) if x.casefold() != 'all'}
types = {x.casefold() for x in filter(None, type_arg.split(','))}
selected = {}
with open(path, newline='') as f:
    reader = csv.DictReader(f)
    required = {'run_number', 'Kin_old', 'target', 'Type', 'prescale'}
    if not required.issubset({x.strip() for x in reader.fieldnames or []}):
        raise SystemExit(f'Config lacks {sorted(required)}')
    for raw in reader:
        row = {k.strip(): (v or '').strip().strip('"') for k, v in raw.items() if k}
        run, kin, target, typ = (row[k] for k in ('run_number', 'Kin_old', 'target', 'Type'))
        if not run.isdigit() or not kin or not target: continue
        if kins and kin not in kins: continue
        if runs and run not in runs: continue
        if targets and target.casefold() not in targets: continue
        if typ.casefold() not in types: continue
        identity = (kin, target, row['prescale'])
        if run in selected and selected[run] != identity:
            raise SystemExit(f'Conflicting metadata for run {run}')
        selected[run] = identity
for run, (kin, target, _) in sorted(selected.items(), key=lambda item: (item[1][0], int(item[0]))):
    if any(c in kin + target for c in '\t\n'):
        raise SystemExit(f'Unsupported tab/newline in run {run}')
    print(run, kin, target, sep='\t')
PY
[[ -s "$work/jobs.tsv" ]] || die "No runs matched selection"
mapfile -t SELECTED_KINS < <(cut -f2 "$work/jobs.tsv" | sort -u)

declare -A SAFE_KIN=() SEEN_SAFE=() RUNS_BY_KIN=() MANIFEST_BY_KEY=() INPUT_BYTES=()
for kin in "${SELECTED_KINS[@]}"; do
  safe="$(safe_name "$kin")"; [[ -n "$safe" ]] || die "Unsafe kinematic: $kin"
  [[ -z "${SEEN_SAFE[$safe]+x}" || "${SEEN_SAFE[$safe]}" == "$kin" ]] || die "Safe-name collision"
  SAFE_KIN[$kin]="$safe"; SEEN_SAFE[$safe]="$kin"
  for prefix in efficiency selection_report; do
    [[ -s "${ANALYSIS_ROOT}/output/efficiency_stuff/${prefix}_${kin}.csv" ]] || die "Missing ${prefix}_${kin}.csv needed for plots"
  done
done

mkdir -p "$TAR_DIR"
source_report="${TAR_DIR}/${WORKFLOW}_run_sources.csv"
[[ ! -e "$source_report" ]] || die "Already exists: $source_report"
source_report_tmp="$(mktemp "${TAR_DIR}/.${WORKFLOW}.XXXXXX.run_sources.csv")"
cleanup_files+=("$source_report_tmp")
printf 'kin,target,run,input_status,source,input_file_count,input_bytes,input_paths\n' > "$source_report_tmp"
shopt -s nullglob
ready=0; missing=0
while IFS=$'\t' read -r run kin target; do
  updated=("${UPDATED_DIR%/}"/nps_hms_coin_"${run}"_*_1_-1.root)
  production=("${PRODUCTION_DIR%/}"/nps_hms_coin_"${run}"_*_1_-1.root)
  if (( ${#updated[@]} )); then source_kind=updated; files=("${updated[@]}")
  elif (( ${#production[@]} )); then source_kind=production; files=("${production[@]}")
  else
    printf '"%s","%s",%s,missing,,0,0,""\n' "$kin" "$target" "$run" >> "$source_report_tmp"
    missing=$((missing + 1)); continue
  fi
  mapfile -t files < <(printf '%s\n' "${files[@]}" | sort -V)
  safe="${SAFE_KIN[$kin]}"; key="${kin}|${run}"
  manifest="${TAR_DIR}/${WORKFLOW}_${safe}_run${run}_inputs.tsv"
  [[ ! -e "$manifest" ]] || die "Already exists: $manifest"
  manifest_tmp="$(mktemp "${TAR_DIR}/.${WORKFLOW}_${safe}_run${run}.XXXXXX.tsv")"
  cleanup_files+=("$manifest_tmp")
  printf 'run\tsource_kind\tinput_access\tbasename\tsource_path\tsize_bytes\n' > "$manifest_tmp"
  paths=""; bytes=0
  for file in "${files[@]}"; do
    [[ -r "$file" ]] || die "Unreadable input: $file"
    if [[ "$INPUT_ACCESS" == stage ]]; then
      perms="$(stat -c %A "$file")"; [[ "${perms:7:1}" == r ]] || die "Staged input is not world-readable: $file"
    fi
    size="$(stat -c %s "$file")"; bytes=$((bytes + size)); basename="$(basename "$file")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$run" "$source_kind" "$INPUT_ACCESS" "$basename" "$file" "$size" >> "$manifest_tmp"
    paths+="${paths:+;}${file}"
  done
  mv "$manifest_tmp" "$manifest"
  # SWIF's staging service reads this file as a service account, not as the
  # submitting login. mktemp creates mode 0600, so publish the finished file.
  chmod 0644 "$manifest"
  printf '"%s","%s",%s,ready,%s,%s,%s,"%s"\n' "$kin" "$target" "$run" "$source_kind" "${#files[@]}" "$bytes" "$paths" >> "$source_report_tmp"
  MANIFEST_BY_KEY[$key]="$manifest"; INPUT_BYTES[$key]="$bytes"; RUNS_BY_KIN[$kin]+="$run "
  ready=$((ready + 1))
done < "$work/jobs.tsv"
chmod 0644 "$source_report_tmp"
mv "$source_report_tmp" "$source_report"
if (( missing )); then
  [[ "$MISSING_INPUT" == skip ]] || die "$missing runs lack replay input; see $source_report"
  echo "[WARN] Skipping $missing runs without replay input; see $source_report" >&2
fi
(( ready )) || die "No usable runs"

active_kins=()
for kin in "${SELECTED_KINS[@]}"; do [[ -n "${RUNS_BY_KIN[$kin]:-}" ]] && active_kins+=("$kin"); done
SELECTED_KINS=("${active_kins[@]}")

# Never mix a new workflow with parts left by an earlier submission.
if (( DRY_RUN == 0 )); then
  for kin in "${SELECTED_KINS[@]}"; do
    safe="${SAFE_KIN[$kin]}"
    [[ ! -e "${OUTPUT_DIR}/${safe}.tar.gz" ]] || die "Result exists: ${OUTPUT_DIR}/${safe}.tar.gz"
    for path in "${OUTPUT_DIR}/parts/${safe}" "${OUTPUT_DIR}/metadata/${safe}/runs"; do
      [[ ! -d "$path" || -z "$(find "$path" -mindepth 1 -maxdepth 1 -print -quit)" ]] || \
        die "Output directory is not empty: $path"
    done
  done
fi

size_bytes() {
  local value="$1" suffix="${1: -1}" number factor=1
  if [[ "$suffix" =~ [kmgt] ]]; then number="${value%?}"; case "$suffix" in k) factor=1024;; m) factor=$((1024**2));; g) factor=$((1024**3));; t) factor=$((1024**4));; esac
  else number="$value"; fi
  echo $((number * factor))
}
if [[ "$INPUT_ACCESS" == stage ]]; then
  allowed="$(size_bytes "$DISK")"
  for key in "${!INPUT_BYTES[@]}"; do (( INPUT_BYTES[$key] * 5 <= allowed * 4 )) || die "$key staged inputs need more --disk (20% headroom required)"; done
fi

# Build once on the submit host; every SWIF job receives this exact executable.
exe="${ANALYSIS_ROOT}/${EXE_REL}"; source_file="${ANALYSIS_ROOT}/${SOURCE_REL}"
if [[ ! -x "$exe" || "$source_file" -nt "$exe" ]] || \
   find "${ANALYSIS_ROOT}/${EFF_REL}" -maxdepth 1 -name '*.h' -newer "$exe" -print -quit | grep -q .; then
  echo "[build] $exe"
  if command -v root-config >/dev/null 2>&1; then
    g++ -std=c++17 -O2 -Wall -Wextra "$source_file" -o "$exe" $(root-config --cflags --libs)
  else
    export NPS_COMPARE_BUILD_SOURCE="$source_file" NPS_COMPARE_BUILD_EXE="$exe"
    csh -f -c 'if ( -f /usr/share/Modules/init/csh ) source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; if ( $status != 0 ) exit 97; g++ -std=c++17 -O2 -Wall -Wextra "$NPS_COMPARE_BUILD_SOURCE" -o "$NPS_COMPARE_BUILD_EXE" `root-config --cflags --libs`'
  fi
fi

runtime_root_name="$(basename "$ANALYSIS_ROOT")"
runtime_tar="${TAR_DIR}/${WORKFLOW}_compare_runtime.tar.gz"; runtime_sha="${runtime_tar}.sha256"
worker_basename="${WORKFLOW}_worker.sh"; worker_script="${TAR_DIR}/${worker_basename}"
[[ ! -e "$runtime_tar" && ! -e "$runtime_sha" && ! -e "$worker_script" ]] || die "Runtime/worker already exists"
archive=("${runtime_root_name}/${EXE_REL}" \
         "${runtime_root_name}/${SOURCE_REL}" "${runtime_root_name}/${PLOT_REL}")
for header in "${ANALYSIS_ROOT}/${EFF_REL}"/*.h; do archive+=("${runtime_root_name}/${EFF_REL}/$(basename "$header")"); done
for kin in "${SELECTED_KINS[@]}"; do
  archive+=("${runtime_root_name}/output/efficiency_stuff/efficiency_${kin}.csv" \
            "${runtime_root_name}/output/efficiency_stuff/selection_report_${kin}.csv")
done
tar_tmp="$(mktemp "${TAR_DIR}/.${WORKFLOW}.XXXXXX.tar.gz")"
cleanup_files+=("$tar_tmp")
tar -C "$(dirname "$ANALYSIS_ROOT")" -czf "$tar_tmp" "${archive[@]}"
tar -tzf "$tar_tmp" >/dev/null
mv "$tar_tmp" "$runtime_tar"
hash="$(sha256sum "$runtime_tar" | awk '{print $1}')"
printf '%s  compare_runtime.tar.gz\n' "$hash" > "$runtime_sha"
worker_tmp="$(mktemp "${TAR_DIR}/.${WORKFLOW}.XXXXXX.worker.sh")"
cleanup_files+=("$worker_tmp")
cp "$SCRIPT_PATH" "$worker_tmp"; chmod 0755 "$worker_tmp"; mv "$worker_tmp" "$worker_script"
chmod 0644 "$runtime_tar" "$runtime_sha"

# Fail before creating workflows if any metadata input cannot be staged.
stage_inputs=("$runtime_tar" "$runtime_sha" "$CONFIG" "$source_report" "$worker_script")
for key in "${!MANIFEST_BY_KEY[@]}"; do stage_inputs+=("${MANIFEST_BY_KEY[$key]}"); done
for file in "${stage_inputs[@]}"; do
  perms="$(stat -c %A "$file")"
  [[ "${perms:7:1}" == r ]] || die "SWIF input is not world-readable: $file ($perms)"
done

if (( DRY_RUN == 0 )); then
  need_cmd swif2
  mkdir -p "$OUTPUT_DIR" "$LOG_DIR"
  for kin in "${SELECTED_KINS[@]}"; do
    safe="${SAFE_KIN[$kin]}"
    mkdir -p "${OUTPUT_DIR}/parts/${safe}" "${OUTPUT_DIR}/metadata/${safe}/runs" "${LOG_DIR}/${safe}"
  done
fi

echo "Workflow base: $WORKFLOW"
echo "Kinematics:   ${SELECTED_KINS[*]}"
echo "Runs:         $ready ready; $missing skipped"
echo "Runtime:      $runtime_tar"
echo "Runtime SHA:  $hash"
echo "Input access: $INPUT_ACCESS"
echo "Output:       $OUTPUT_DIR"
echo "Logs:         $LOG_DIR"

swif() { if (( DRY_RUN )); then printf '[dry-run]'; printf ' %q' swif2 "$@"; printf '\n'; else swif2 "$@"; fi; }
workflows=()
for kin in "${SELECTED_KINS[@]}"; do
  safe="${SAFE_KIN[$kin]}"; workflow="${WORKFLOW}_${safe}"; workflows+=("$workflow")
  read -r -a kin_runs <<< "${RUNS_BY_KIN[$kin]}"
  swif create "$workflow" -max-concurrent "$MAX_CONCURRENT"
  antecedents=()
  for run in "${kin_runs[@]}"; do
    key="${kin}|${run}"; job="run_${run}"; antecedents+=(-antecedent "$job")
    inputs=(-input compare_runtime.tar.gz "file:${runtime_tar}" \
            -input compare_runtime.tar.gz.sha256 "file:${runtime_sha}" \
            -input compare_config.csv "file:${CONFIG}" \
            -input input_manifest.tsv "file:${MANIFEST_BY_KEY[$key]}" \
            -input run_source_report.csv "file:${source_report}" \
            -input "$worker_basename" "file:${worker_script}")
    if [[ "$INPUT_ACCESS" == stage ]]; then
      while IFS=$'\t' read -r row_run source_kind access basename source_path bytes; do
        [[ "$row_run" == run ]] && continue
        inputs+=(-input "$basename" "file:${source_path}")
      done < "${MANIFEST_BY_KEY[$key]}"
    fi
    swif add-job "$workflow" -name "$job" -account "$ACCOUNT" -partition "$PARTITION" \
      -cores "$CORES" -ram "$RAM" -disk "$DISK" -time "$WALLTIME" \
      -stdout "${LOG_DIR}/${safe}/${job}.out" -stderr "${LOG_DIR}/${safe}/${job}.err" \
      "${inputs[@]}" \
      -output compare_livetimes_part.csv "file:${OUTPUT_DIR}/parts/${safe}/run${run}.csv" \
      -output compare_livetimes_job_manifest.txt "file:${OUTPUT_DIR}/metadata/${safe}/runs/run${run}.txt" \
      /bin/bash "$worker_basename" --run-one compare_runtime.tar.gz "$runtime_root_name" "$kin" "$run"
  done
  result="${OUTPUT_DIR}/${safe}.tar.gz"
  swif add-job "$workflow" -name "finalize_${safe}" "${antecedents[@]}" \
    -account "$ACCOUNT" -partition "$PARTITION" -cores 1 -ram "$RAM" -disk "$DISK" -time "$WALLTIME" \
    -stdout "${LOG_DIR}/${safe}/finalize.out" -stderr "${LOG_DIR}/${safe}/finalize.err" \
    -input compare_runtime.tar.gz "file:${runtime_tar}" \
    -input compare_runtime.tar.gz.sha256 "file:${runtime_sha}" \
    -input compare_config.csv "file:${CONFIG}" \
    -input run_source_report.csv "file:${source_report}" \
    -input "$worker_basename" "file:${worker_script}" \
    /bin/bash "$worker_basename" --finalize-one compare_runtime.tar.gz "$runtime_root_name" "$kin" \
      "$OUTPUT_DIR" "$result" "${OUTPUT_DIR}/parts/${safe}"
done

if (( CREATE_ONLY == 0 )); then for workflow in "${workflows[@]}"; do swif run "$workflow"; done; fi
if (( DRY_RUN )); then echo "Dry run complete; SWIF2 state unchanged."
elif (( CREATE_ONLY )); then printf 'Created, not started: %s\n' "${workflows[@]}"
else printf 'Submitted: swif2 status %s\n' "${workflows[@]}"; fi
