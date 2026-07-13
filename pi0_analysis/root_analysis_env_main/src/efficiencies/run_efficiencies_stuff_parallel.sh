#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

EXE="${SCRIPT_DIR}/compute_efficiencies_stuff"
CONFIG_CSV="/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/config/nps_dvcs_all_kins_main.csv"
UPDATED_DIR="/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/updated"
PRODUCTION_DIR="/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/production"
OUTPUT_DIR="/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/output/efficiency_stuff"
TYPES_CSV="production,Production"
JOBS="$(nproc)"

ALL_KINS=0
KIN_LIST=()
RUN_FILTER=""

print_help() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --exe <path>              Path to compute_efficiencies_stuff executable
  --config <path>           Config CSV path
  --updated-dir <path>      Updated ROOT directory
  --production-dir <path>   Production ROOT directory
  --output-dir <path>       Final output directory
  --types <a,b,c>           Allowed Type values (default: production,Production)
  --kin <Kin_old>           Restrict to one kinematic setting (repeatable)
  --all-kins                Process all kinematic settings
  --run <run_number>        Restrict to one run number
  --jobs <N>                Parallel workers (default: nproc)
  --help                    Show this message
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --exe)
      EXE="$2"
      shift 2
      ;;
    --config)
      CONFIG_CSV="$2"
      shift 2
      ;;
    --updated-dir)
      UPDATED_DIR="$2"
      shift 2
      ;;
    --production-dir)
      PRODUCTION_DIR="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --types)
      TYPES_CSV="$2"
      shift 2
      ;;
    --kin)
      KIN_LIST+=("$2")
      shift 2
      ;;
    --all-kins)
      ALL_KINS=1
      shift
      ;;
    --run)
      RUN_FILTER="$2"
      shift 2
      ;;
    --jobs)
      JOBS="$2"
      shift 2
      ;;
    --help|-h)
      print_help
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      print_help
      exit 1
      ;;
  esac
done

if [[ ${ALL_KINS} -eq 1 && ${#KIN_LIST[@]} -gt 0 ]]; then
  echo "Use either --all-kins or --kin, not both." >&2
  exit 1
fi

if [[ ! -x "${EXE}" ]]; then
  if command -v root-config >/dev/null 2>&1; then
    SRC="${SCRIPT_DIR}/compute_efficiencies_stuff.cxx"
    echo "[build] Building executable: ${EXE}"
    g++ -O3 -march=native -std=c++17 "${SRC}" -o "${EXE}" $(root-config --cflags --libs)
  else
    echo "Executable not found and root-config is unavailable: ${EXE}" >&2
    exit 1
  fi
fi

if [[ ! -f "${CONFIG_CSV}" ]]; then
  echo "Config CSV not found: ${CONFIG_CSV}" >&2
  exit 1
fi

sanitize_name() {
  local input="$1"
  echo "${input}" | sed 's/[^[:alnum:]]/_/g'
}

join_by_pipe() {
  local IFS='|'
  echo "$*"
}

render_progress_bar() {
  local done="$1"
  local total="$2"
  local width="${3:-40}"

  if [[ "${total}" -le 0 ]]; then
    printf "[----------------------------------------]   0%% (0/0)"
    return
  fi

  if [[ "${done}" -gt "${total}" ]]; then
    done="${total}"
  fi

  local percent=$(( done * 100 / total ))
  local filled=$(( done * width / total ))
  local empty=$(( width - filled ))

  local filled_bar empty_bar
  filled_bar="$(printf "%*s" "${filled}" "" | tr ' ' '#')"
  empty_bar="$(printf "%*s" "${empty}" "" | tr ' ' '-')"

  printf "[%s%s] %3d%% (%d/%d)" "${filled_bar}" "${empty_bar}" "${percent}" "${done}" "${total}"
}

trim_ws() {
  local s="$1"
  s="${s#${s%%[![:space:]]*}}"
  s="${s%${s##*[![:space:]]}}"
  echo "$s"
}

get_all_kins() {
  awk -F',' '
    function trim(s) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
      return s
    }
    NR==1 {
      for (i=1; i<=NF; ++i) {
        if (trim($i) == "Kin_old") kin_col=i
      }
      next
    }
    NR>1 && kin_col > 0 {
      kin = trim($kin_col)
      if (kin != "") seen[kin]=1
    }
    END {
      for (k in seen) print k
    }' "${CONFIG_CSV}" | sort
}

prompt_for_kinematics() {
  local -a all_kins=("$@")
  if [[ ${#all_kins[@]} -eq 0 ]]; then
    echo "No kinematic settings found in CSV: ${CONFIG_CSV}" >&2
    return 1
  fi

  while true; do
    echo "Available kinematic settings (Kin_old):"
    local i
    for ((i=0; i<${#all_kins[@]}; ++i)); do
      printf "  %3d) %s\n" "$((i + 1))" "${all_kins[$i]}"
    done
    echo ""
    echo "Choose one of the following:"
    echo "  - all"
    echo "  - one or more numbers separated by commas (example: 1,4,7)"
    echo "  - one or more Kin_old names separated by commas"
    read -r -p "Selection [all]: " user_input
    user_input="$(trim_ws "${user_input}")"
    [[ -z "${user_input}" ]] && user_input="all"

    if [[ "${user_input}" == "all" ]]; then
      SELECTED_KINS=("${all_kins[@]}")
      echo "[select] Using all ${#SELECTED_KINS[@]} kinematic settings."
      return 0
    fi

    local -a tokens=()
    IFS=',' read -r -a tokens <<< "${user_input}"

    local -a selected=()
    local token matched idx kin
    declare -A seen=()

    for token in "${tokens[@]}"; do
      token="$(trim_ws "${token}")"
      [[ -z "${token}" ]] && continue

      matched=0
      if [[ "${token}" =~ ^[0-9]+$ ]]; then
        idx=$((token - 1))
        if (( idx >= 0 && idx < ${#all_kins[@]} )); then
          kin="${all_kins[$idx]}"
          matched=1
        fi
      else
        for kin in "${all_kins[@]}"; do
          if [[ "${kin}" == "${token}" ]]; then
            matched=1
            break
          fi
        done
      fi

      if (( matched == 1 )); then
        if [[ -z "${seen[${kin}]+x}" ]]; then
          selected+=("${kin}")
          seen["${kin}"]=1
        fi
      else
        echo "[warn] Invalid selection token: ${token}"
      fi
    done

    if [[ ${#selected[@]} -gt 0 ]]; then
      SELECTED_KINS=("${selected[@]}")
      echo "[select] Chosen kinematic settings: ${SELECTED_KINS[*]}"
      return 0
    fi

    echo "[select] No valid kinematic choice provided. Please try again."
    echo ""
  done
}

SELECTED_KINS=()
if [[ ${ALL_KINS} -eq 1 ]]; then
  mapfile -t SELECTED_KINS < <(get_all_kins)
elif [[ ${#KIN_LIST[@]} -gt 0 ]]; then
  SELECTED_KINS=("${KIN_LIST[@]}")
else
  mapfile -t AVAILABLE_KINS < <(get_all_kins)
  if [[ ! -t 0 ]]; then
    echo "[select] No --kin/--all-kins provided and no interactive terminal detected; defaulting to all kinematics."
    SELECTED_KINS=("${AVAILABLE_KINS[@]}")
  else
    prompt_for_kinematics "${AVAILABLE_KINS[@]}"
  fi
fi

if [[ ${#SELECTED_KINS[@]} -eq 0 ]]; then
  echo "No kinematic settings selected." >&2
  exit 1
fi

KIN_PATTERN="$(join_by_pipe "${SELECTED_KINS[@]}")"
TYPE_ITEMS=()
IFS=',' read -r -a TYPE_ITEMS <<< "${TYPES_CSV}"
TYPE_PATTERN="$(join_by_pipe "${TYPE_ITEMS[@]}")"

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

mkdir -p "${OUTPUT_DIR}"
RUN_LOG_DIR="${OUTPUT_DIR}/run_logs"
mkdir -p "${RUN_LOG_DIR}"

JOB_LIST="${TMP_DIR}/jobs.txt"

awk -F',' -v kin_pattern="${KIN_PATTERN}" -v type_pattern="${TYPE_PATTERN}" -v run_filter="${RUN_FILTER}" '
  function trim(s) {
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
    return s
  }
  BEGIN {
    n_kin = split(kin_pattern, kin_parts, "|")
    for (i=1; i<=n_kin; ++i) {
      kin_parts[i] = trim(kin_parts[i])
      if (kin_parts[i] != "") wanted_kin[kin_parts[i]] = 1
    }

    n_type = split(type_pattern, type_parts, "|")
    for (i=1; i<=n_type; ++i) {
      type_parts[i] = trim(type_parts[i])
      if (type_parts[i] != "") wanted_type[type_parts[i]] = 1
    }
  }
  NR == 1 {
    for (i=1; i<=NF; ++i) {
      name = trim($i)
      if (name == "run_number") run_col=i
      if (name == "Kin_old") kin_col=i
      if (name == "Type") type_col=i
    }
    next
  }
  NR > 1 {
    if (!(kin_col > 0 && run_col > 0 && type_col > 0)) next
    kin = trim($kin_col)
    run = trim($run_col)
    type = trim($type_col)

    if (!(kin in wanted_kin)) next
    if (!(type in wanted_type)) next
    if (run !~ /^[0-9]+$/) next
    if (run_filter != "" && run != run_filter) next

    key = kin "," run
    if (!(key in seen)) {
      print key
      seen[key] = 1
    }
  }' "${CONFIG_CSV}" | sort > "${JOB_LIST}"

TOTAL_JOBS="$(wc -l < "${JOB_LIST}")"
if [[ "${TOTAL_JOBS}" -eq 0 ]]; then
  echo "No runs matched selected kinematics/types." >&2
  exit 1
fi

if ! [[ "${JOBS}" =~ ^[0-9]+$ ]] || [[ "${JOBS}" -lt 1 ]]; then
  echo "[warn] Invalid --jobs value '${JOBS}', defaulting to 1"
  JOBS=1
fi

EFFECTIVE_JOBS="${JOBS}"
if [[ "${EFFECTIVE_JOBS}" -gt "${TOTAL_JOBS}" ]]; then
  EFFECTIVE_JOBS="${TOTAL_JOBS}"
fi

echo "[parallel] Launching ${TOTAL_JOBS} run jobs with ${EFFECTIVE_JOBS} workers"

PROGRESS_FILE="${TMP_DIR}/progress.done"
START_FILE="${TMP_DIR}/progress.start"
: > "${PROGRESS_FILE}"
: > "${START_FILE}"

export EXE CONFIG_CSV UPDATED_DIR PRODUCTION_DIR TYPES_CSV TMP_DIR PROGRESS_FILE START_FILE RUN_LOG_DIR

xargs -P "${EFFECTIVE_JOBS}" -I {} bash -c '
  line="$1"
  IFS="," read -r kin run <<< "${line}"

  safe_kin="$(echo "${kin}" | sed "s/[^[:alnum:]]/_/g")"
  job_out="${TMP_DIR}/jobs/${safe_kin}/${run}"
  mkdir -p "${job_out}"
  run_log="${RUN_LOG_DIR}/${safe_kin}_run${run}.log"
  run_debug_log="${RUN_LOG_DIR}/${safe_kin}_run${run}_debug.log"

  printf "%s,%s\n" "${kin}" "${run}" >> "${START_FILE}"

  status=0
  "${EXE}" \
    --config "${CONFIG_CSV}" \
    --updated-dir "${UPDATED_DIR}" \
    --production-dir "${PRODUCTION_DIR}" \
    --output-dir "${job_out}" \
    --kin "${kin}" \
    --types "${TYPES_CSV}" \
    --run "${run}" \
    --log "${run_debug_log}" \
    --no-interactive \
    > "${run_log}" 2>&1 || status=$?

  printf "%s,%s,%s\n" "${kin}" "${run}" "${status}" >> "${PROGRESS_FILE}"
  exit "${status}"
' _ {} < "${JOB_LIST}" &

XARGS_PID=$!
XARGS_STATUS=0

if [[ -t 1 ]]; then
  SPIN_CHARS='|/-\\'
  SPIN_INDEX=0
  while kill -0 "${XARGS_PID}" 2>/dev/null; do
    STARTED_JOBS="$(wc -l < "${START_FILE}")"
    DONE_JOBS="$(wc -l < "${PROGRESS_FILE}")"
    RUNNING_JOBS=$(( STARTED_JOBS - DONE_JOBS ))
    if [[ "${RUNNING_JOBS}" -lt 0 ]]; then
      RUNNING_JOBS=0
    fi

    SPIN_CHAR="${SPIN_CHARS:${SPIN_INDEX}:1}"
    SPIN_INDEX=$(( (SPIN_INDEX + 1) % 4 ))

    printf "\r[progress] "
    render_progress_bar "${DONE_JOBS}" "${TOTAL_JOBS}" 40
    printf "  %s started:%d running:%d done:%d/%d" \
      "${SPIN_CHAR}" "${STARTED_JOBS}" "${RUNNING_JOBS}" "${DONE_JOBS}" "${TOTAL_JOBS}"
    sleep 1
  done
fi

wait "${XARGS_PID}" || XARGS_STATUS=$?

STARTED_JOBS="$(wc -l < "${START_FILE}")"
DONE_JOBS="$(wc -l < "${PROGRESS_FILE}")"
RUNNING_JOBS=$(( STARTED_JOBS - DONE_JOBS ))
if [[ "${RUNNING_JOBS}" -lt 0 ]]; then
  RUNNING_JOBS=0
fi
if [[ -t 1 ]]; then
  printf "\r[progress] "
  render_progress_bar "${DONE_JOBS}" "${TOTAL_JOBS}" 40
  printf "  done started:%d running:%d done:%d/%d" \
    "${STARTED_JOBS}" "${RUNNING_JOBS}" "${DONE_JOBS}" "${TOTAL_JOBS}"
  printf "\n"
else
  echo "[progress] started ${STARTED_JOBS}, completed ${DONE_JOBS}/${TOTAL_JOBS} jobs"
fi

if [[ "${XARGS_STATUS}" -ne 0 ]]; then
  FAILED_JOBS="$(awk -F',' '$3 != 0 { c++ } END { print (c+0) }' "${PROGRESS_FILE}")"
  echo "[error] ${FAILED_JOBS} run job(s) failed. Check per-run logs under ${RUN_LOG_DIR}." >&2
  exit "${XARGS_STATUS}"
fi

for kin in "${SELECTED_KINS[@]}"; do
  safe_kin="$(sanitize_name "${kin}")"
  final_csv="${OUTPUT_DIR}/efficiency_${safe_kin}.csv"
  final_selection_csv="${OUTPUT_DIR}/selection_report_${safe_kin}.csv"

  mapfile -t parts < <(find "${TMP_DIR}/jobs/${safe_kin}" -type f -name "efficiency_${safe_kin}.csv" | sort)
  if [[ ${#parts[@]} -eq 0 ]]; then
    echo "[merge] No partial CSV files found for ${kin}; skipping"
    continue
  fi

  head -n 1 "${parts[0]}" > "${final_csv}"
  for part in "${parts[@]}"; do
    tail -n +2 "${part}" >> "${final_csv}"
  done

  {
    head -n 1 "${final_csv}"
    tail -n +2 "${final_csv}" | sort -t',' -k1,1n
  } > "${final_csv}.tmp"
  mv "${final_csv}.tmp" "${final_csv}"

  echo "[merge] Wrote ${final_csv}"

  mapfile -t selection_parts < <(find "${TMP_DIR}/jobs/${safe_kin}" -type f -name "selection_report_${safe_kin}.csv" | sort)
  if [[ ${#selection_parts[@]} -gt 0 ]]; then
    head -n 1 "${selection_parts[0]}" > "${final_selection_csv}"
    for part in "${selection_parts[@]}"; do
      tail -n +2 "${part}" >> "${final_selection_csv}"
    done

    {
      head -n 1 "${final_selection_csv}"
      tail -n +2 "${final_selection_csv}" | sort -t',' -k1,1n -k4,4n
    } > "${final_selection_csv}.tmp"
    mv "${final_selection_csv}.tmp" "${final_selection_csv}"

    echo "[merge] Wrote ${final_selection_csv}"
  fi
done

SUMMARY_CSV="${OUTPUT_DIR}/efficiency_runs_processed.csv"
echo "kinematic_setting,total_runs_listed,production_runs,runs_selected_by_type,runs_processed,runs_not_found,malformed_rows_skipped" > "${SUMMARY_CSV}"

for kin in "${SELECTED_KINS[@]}"; do
  safe_kin="$(sanitize_name "${kin}")"
  merged_csv="${OUTPUT_DIR}/efficiency_${safe_kin}.csv"

	  total_runs_listed="$(awk -F',' -v kin="${kin}" '
      function trim(s) {
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
        return s
      }
	    NR==1 {
	      for (i=1; i<=NF; ++i) {
	        name = trim($i)
	        if (name == "run_number") run_col=i
	        if (name == "Kin_old") kin_col=i
	      }
	      next
	    }
	    NR>1 && kin_col>0 && run_col>0 {
	      row_kin = trim($kin_col)
	      run = trim($run_col)
	      if (row_kin == kin && run ~ /^[0-9]+$/) seen[run]=1
	    }
    END {
      c=0
      for (x in seen) c++
      print c
    }' "${CONFIG_CSV}")"

	  production_runs="$(awk -F',' -v kin="${kin}" '
      function trim(s) {
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
        return s
      }
	    NR==1 {
	      for (i=1; i<=NF; ++i) {
	        name = trim($i)
	        if (name == "run_number") run_col=i
	        if (name == "Kin_old") kin_col=i
	        if (name == "Type") type_col=i
	      }
	      next
	    }
	    NR>1 && kin_col>0 && run_col>0 && type_col>0 {
	      row_kin = trim($kin_col)
	      type = trim($type_col)
	      run = trim($run_col)
	      if (row_kin == kin && (type == "production" || type == "Production") && run ~ /^[0-9]+$/) seen[run]=1
	    }
    END {
      c=0
      for (x in seen) c++
      print c
    }' "${CONFIG_CSV}")"

	  runs_selected_by_type="$(awk -F',' -v kin="${kin}" -v type_pattern="${TYPE_PATTERN}" -v run_filter="${RUN_FILTER}" '
      function trim(s) {
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
        return s
      }
	    BEGIN {
	      n_type = split(type_pattern, type_parts, "|")
	      for (i=1; i<=n_type; ++i) {
	        type_parts[i] = trim(type_parts[i])
	        if (type_parts[i] != "") wanted[type_parts[i]]=1
	      }
	    }
	    NR==1 {
	      for (i=1; i<=NF; ++i) {
	        name = trim($i)
	        if (name == "run_number") run_col=i
	        if (name == "Kin_old") kin_col=i
	        if (name == "Type") type_col=i
	      }
	      next
	    }
	    NR>1 && kin_col>0 && run_col>0 && type_col>0 {
	      row_kin = trim($kin_col)
	      type = trim($type_col)
	      run = trim($run_col)
	      if (row_kin != kin) next
	      if (!(type in wanted)) next
	      if (run !~ /^[0-9]+$/) next
	      if (run_filter != "" && run != run_filter) next
	      seen[run]=1
	    }
    END {
      c=0
      for (x in seen) c++
      print c
    }' "${CONFIG_CSV}")"

  runs_processed="0"
  runs_not_found="0"
  if [[ -f "${merged_csv}" ]]; then
    status_col="$(awk -F',' 'NR==1 {
      for (i=1; i<=NF; ++i) {
        name=$i
        gsub(/"/, "", name)
        if (name == "run_processing_status") {
          print i
          exit
        }
      }
    }' "${merged_csv}")"

    if [[ -n "${status_col}" ]]; then
      runs_processed="$(awk -F',' -v sc="${status_col}" 'NR>1 {
        status=$sc
        gsub(/"/, "", status)
        if (status == "processed" || status == "processed_partial") c++
      } END { print (c+0) }' "${merged_csv}")"

      runs_not_found="$(awk -F',' -v sc="${status_col}" 'NR>1 {
        status=$sc
        gsub(/"/, "", status)
        if (status == "missing_root_files") c++
      } END { print (c+0) }' "${merged_csv}")"
    fi
  fi

  echo "\"${kin}\",${total_runs_listed},${production_runs},${runs_selected_by_type},${runs_processed},${runs_not_found},0" >> "${SUMMARY_CSV}"
done

echo "[done] Summary written to ${SUMMARY_CSV}"
