#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

ROOT_MACRO="${SCRIPT_DIR}/nps_analysis_main.C"
ACCEPTANCE_CUTS_IMPL="${SCRIPT_DIR}/acceptance_cuts.cpp"
COMBINE_SCRIPT="${SCRIPT_DIR}/combine_analysis_branches.py"
CONFIG_CSV="${REPO_ROOT}/config/nps_dvcs_all_kins_main.csv"
OUTPUT_BASE="${REPO_ROOT}/output"
ACCEPTANCE_CUTS_CONFIG="${NPS_ACCEPTANCE_CUTS_CONFIG:-${REPO_ROOT}/config/acceptance_cuts.conf}"

UPDATED_DIR="/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/updated"
PRODUCTION_DIR="/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/production"
WAVEFORM_DIR="/lustre24/expphy/volatile/hallc/nps/hhuang/farmFile/Production/DVCS"

ROOT_CMD="${ROOT_CMD:-root}"
PYTHON_CMD="${PYTHON_CMD:-python3}"
MODE="auto"
SOURCE=""
INPUT_DIR=""
TYPES_CSV="production,Production"
JOBS="$(nproc)"
TIMEOUT_SEC="0"
RUN_FILTER=""
RUN_FILTERS=()
USE_GEVNUM_CUT="${NPS_USE_GEVNUM_CUT:-ask}"
COMBINE_AFTER_RUN="yes"
COMBINE_TARGET="LH2"

RUN_SMEARING_STAGE="no"
SMEAR_SCRIPT="${REPO_ROOT}/src/simulation_smearing/run_smearing_pipeline.sh"
SMEAR_TARGET=""
SMEAR_MODE_HINT=""
SMEAR_ROOT_DIR=""
SMEAR_COMBINED_FILE=""
SMEAR_SIM_FILE=""
SMEAR_SIM_SMEARED_FILE=""
SMEAR_OUT_FILE=""
SMEAR_SECTION_MAP_FILE=""
SMEAR_INTERP_FILE=""
SMEAR_DATA_TREE=""
SMEAR_SIM_TREE=""
SMEAR_NX=""
SMEAR_NY=""
SMEAR_X_MIN=""
SMEAR_X_MAX=""
SMEAR_Y_MIN=""
SMEAR_Y_MAX=""
SMEAR_OVERLAP=""
SMEAR_NSMEAR=""
SMEAR_SIM_EXCL_INPUT=""
SMEAR_SIM_SIDIS_INPUT=""
SMEAR_SIM_DELTA_INPUT=""

RUN_XSEC_STAGE="no"
XSEC_SCRIPT="${REPO_ROOT}/src/xsec_extract/run_xsec_pipeline.sh"
XSEC_TARGET=""
XSEC_ROOT_DIR=""
XSEC_DATA_FILE=""
XSEC_SIM_FILE=""
XSEC_OUT_DIR=""
XSEC_OUT_ROOT=""
XSEC_OUT_CSV=""
XSEC_OUT_SLICE_CSV=""
XSEC_ALL_PLOTS_PDF=""

ALL_KINS=0
KIN_LIST=()
SMEAR_KINS=()
XSEC_KINS=()

CSV_HEADER="run,accumulated_charge(mC),hel_pos_charge(mC),hel_neg_charge(mC),current_mean_uA,Beam_Time(s),total_entries,pass_hms,pass_hms_nps,total_coin_entries,estimated_time_accidentals,chi2_ndf_comb_bg,pi0_mu_MeV,pi0_sigma_MeV,pi0_signal_counts,mmiss_p_mean_GeV,mmiss_p_sigma_GeV,s1x_peak,s1x_err,s1y_peak,s1y_err,s2x_peak,s2x_err,s2y_peak,s2y_err,run_status"

print_help() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --config <path>           Main kinematic config CSV
  --acceptance-config <path>
                            Acceptance-cuts config used by analysis/smearing stages
  --output-base <path>      Canonical output base directory
  --updated-dir <path>      Updated HCANA ROOT directory
  --production-dir <path>   Production HCANA ROOT directory
  --waveform-dir <path>     Waveform ROOT directory
  --input-dir <path>        Explicit input ROOT directory (overrides source dirs)
  --source <updated|production|waveform>
                            Force one source directory (otherwise auto-probes per run)
  --mode <auto|hcana|waveform>
                            Analysis mode (default: auto)
  --types <a,b,c>           Allowed Type values from config CSV
  --no-combine              Skip post-processing combine step
  --combine-target <name>   Target to combine (default: LH2)
  --run-smearing            Run simulation smearing stage after combine
  --smearing-kin <Kin_old>  Restrict smearing to one Kin_old (repeatable)
  --smearing-script <path>  Smearing pipeline script path
  --smearing-target <name>  Smearing target token (default: combine target)
  --smearing-mode <mode>    Mode hint passed to smearing SIM producer (default: analysis --mode)
  --smearing-root-dir <path>            Override smearing root directory
  --smearing-combined-file <path>       Override smearing combined ROOT path
  --smearing-sim-file <path>            Override unsmeared SIM ROOT path
  --smearing-sim-smeared-file <path>    Override smeared SIM ROOT path
  --smearing-out-file <path>            Override smearing fit ROOT output
  --smearing-section-map-file <path>    Override section_map.csv output path
  --smearing-interp-file <path>         Override interpolated ROOT output path
  --smearing-data-tree <name>           Override data tree name for smearing fit
  --smearing-sim-tree <name>            Override simulation tree name for smearing fit
  --smearing-nx <int>                   Override smearing grid nx
  --smearing-ny <int>                   Override smearing grid ny
  --smearing-x-min <float>              Override smearing x-min
  --smearing-x-max <float>              Override smearing x-max
  --smearing-y-min <float>              Override smearing y-min
  --smearing-y-max <float>              Override smearing y-max
  --smearing-overlap <float>            Override section overlap fraction
  --smearing-nsmear <int>               Override smearing iterations/event
  --smearing-sim-excl-input <path>      Override exclusive SIM input for producer
  --smearing-sim-sidis-input <path>     Override SIDIS SIM input for producer
  --smearing-sim-delta-input <path>     Override delta pi0 SIM input for producer
  --run-xsec                Run xsec extraction stage after combine/smearing
  --xsec-kin <Kin_old>      Restrict xsec stage to one Kin_old (repeatable)
  --xsec-script <path>      Xsec pipeline script path
  --xsec-target <name>      Xsec target token (default: combine target)
  --xsec-root-dir <path>    Override xsec root directory
  --xsec-data-file <path>   Override xsec data ROOT input path
  --xsec-sim-file <path>    Override xsec simulation ROOT input path
  --xsec-out-dir <path>     Override xsec output directory
  --xsec-out-root <path>    Override xsec output ROOT path
  --xsec-out-csv <path>     Override xsec output summary CSV path
  --xsec-out-slice-csv <path>
                            Override xsec output slice CSV path
  --xsec-all-plots-pdf <path>
                            Override xsec combined plots PDF path
  --kin <Kin_old>           Restrict to one kinematic setting (repeatable)
  --all-kins                Process all kinematic settings
  --run <run_number...>     Restrict to one or more runs (space- or comma-separated)
  --gevnum-cut <ask|yes|no> Use optional g.evnum event cut (default: ask)
  --jobs <N>                Parallel workers (default: nproc)
  --timeout <seconds>       Per-run timeout (0 disables timeout)
  --help                    Show this message
EOF
}

trim_ws() {
  local s="$1"
  s="${s#${s%%[![:space:]]*}}"
  s="${s%${s##*[![:space:]]}}"
  echo "$s"
}

sanitize_name() {
  echo "$1" | sed 's/[^[:alnum:]_-]/_/g'
}

combine_cut_debug_pdfs() {
  local kin="$1"
  local safe_kin plots_dir out_pdf
  safe_kin="$(sanitize_name "${kin}")"
  plots_dir="${OUTPUT_BASE}/${safe_kin}/plots"
  out_pdf="${plots_dir}/cut_debug_${safe_kin}.pdf"

  [[ -d "${plots_dir}" ]] || return 0

  local -a pdfs=()
  mapfile -t pdfs < <(find "${plots_dir}" -maxdepth 1 -type f -name 'cut_debug_run*.pdf' | sort -V)
  if [[ ${#pdfs[@]} -eq 0 ]]; then
    echo "[cut-debug] No per-run cut-debug PDFs found for kin=${kin}"
    return 0
  fi

  if command -v pdfunite >/dev/null 2>&1; then
    pdfunite "${pdfs[@]}" "${out_pdf}"
    echo "[cut-debug] Wrote ${out_pdf} from ${#pdfs[@]} run PDF(s)"
  elif command -v gs >/dev/null 2>&1; then
    gs -q -dBATCH -dNOPAUSE -sDEVICE=pdfwrite -sOutputFile="${out_pdf}" "${pdfs[@]}"
    echo "[cut-debug] Wrote ${out_pdf} from ${#pdfs[@]} run PDF(s)"
  else
    echo "[cut-debug] Could not merge cut-debug PDFs for kin=${kin}: install pdfunite or ghostscript. Per-run PDFs remain in ${plots_dir}." >&2
  fi
}

join_by_pipe() {
  local IFS='|'
  echo "$*"
}

get_all_kins() {
  awk -F',' '
    function trim(s) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
      return s
    }
    NR==1 {
      for (i=1; i<=NF; ++i) {
        h = trim($i)
        if (h == "Kin_old") kin_col=i
      }
      next
    }
    NR>1 && kin_col > 0 {
      v = trim($kin_col)
      if (v != "") seen[v]=1
    }
    END {
      for (k in seen) print k
    }' "${CONFIG_CSV}" | sort
}

get_kins_for_run() {
  local run="$1"
  awk -F',' -v run="${run}" '
    function trim(s) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
      return s
    }
    NR == 1 {
      for (i=1; i<=NF; ++i) {
        h = trim($i)
        if (h == "run_number") run_col=i
        if (h == "Kin_old") kin_col=i
      }
      next
    }
    NR > 1 {
      if (!(run_col > 0 && kin_col > 0)) next
      rv = trim($run_col)
      kv = trim($kin_col)
      if (rv == run && kv != "") seen[kv] = 1
    }
    END {
      for (k in seen) print k
    }' "${CONFIG_CSV}" | sort
}

prompt_for_kinematics() {
  local -a all_kins=("$@")
  if [[ ${#all_kins[@]} -eq 0 ]]; then
    echo "No kinematic settings found in CSV." >&2
    return 1
  fi

  while true; do
    echo "Available kinematic settings (Kin_old):"
    local i
    for ((i=0; i<${#all_kins[@]}; ++i)); do
      printf "  %3d) %s\n" "$((i + 1))" "${all_kins[$i]}"
    done
    echo
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
    echo
  done
}

prompt_for_source() {
  while true; do
    echo "Select input ROOT source:"
    echo "  1) updated   (HCANA files from updated replay directory)"
    echo "  2) production (HCANA files from production replay directory)"
    echo "  3) waveform  (waveform-analyzed ROOT files)"
    read -r -p "Selection [1]: " source_choice
    source_choice="$(trim_ws "${source_choice}")"
    [[ -z "${source_choice}" ]] && source_choice="1"

    case "${source_choice}" in
      1|updated)
        SOURCE="updated"
        INPUT_DIR="${UPDATED_DIR}"
        [[ "${MODE}" == "auto" ]] && MODE="hcana"
        return 0
        ;;
      2|production)
        SOURCE="production"
        INPUT_DIR="${PRODUCTION_DIR}"
        [[ "${MODE}" == "auto" ]] && MODE="hcana"
        return 0
        ;;
      3|waveform)
        SOURCE="waveform"
        INPUT_DIR="${WAVEFORM_DIR}"
        [[ "${MODE}" == "auto" ]] && MODE="waveform"
        return 0
        ;;
      *)
        echo "[warn] Invalid source choice: ${source_choice}"
        ;;
    esac
  done
}

prompt_for_gevnum_cut() {
  while true; do
    read -r -p "Apply optional g.evnum good-event cut from selection report CSV? [Y/n]: " cut_choice
    cut_choice="$(trim_ws "${cut_choice}")"
    [[ -z "${cut_choice}" ]] && cut_choice="y"
    cut_choice="$(echo "${cut_choice}" | tr '[:upper:]' '[:lower:]')"

    case "${cut_choice}" in
      y|yes)
        USE_GEVNUM_CUT="yes"
        return 0
        ;;
      n|no)
        USE_GEVNUM_CUT="no"
        return 0
        ;;
      *)
        echo "[warn] Invalid choice: ${cut_choice}"
        ;;
    esac
  done
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_CSV="$2"
      shift 2
      ;;
    --acceptance-config)
      ACCEPTANCE_CUTS_CONFIG="$2"
      shift 2
      ;;
    --output-base)
      OUTPUT_BASE="$2"
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
    --waveform-dir)
      WAVEFORM_DIR="$2"
      shift 2
      ;;
    --input-dir)
      INPUT_DIR="$2"
      shift 2
      ;;
    --source)
      SOURCE="$2"
      shift 2
      ;;
    --mode)
      MODE="$2"
      shift 2
      ;;
    --types)
      TYPES_CSV="$2"
      shift 2
      ;;
    --no-combine)
      COMBINE_AFTER_RUN="no"
      shift
      ;;
    --combine-target)
      COMBINE_TARGET="$2"
      shift 2
      ;;
    --run-smearing)
      RUN_SMEARING_STAGE="yes"
      shift
      ;;
    --smearing-kin)
      SMEAR_KINS+=("$2")
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-script)
      SMEAR_SCRIPT="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-target)
      SMEAR_TARGET="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-mode)
      SMEAR_MODE_HINT="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-root-dir)
      SMEAR_ROOT_DIR="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-combined-file)
      SMEAR_COMBINED_FILE="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-sim-file)
      SMEAR_SIM_FILE="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-sim-smeared-file)
      SMEAR_SIM_SMEARED_FILE="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-out-file)
      SMEAR_OUT_FILE="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-section-map-file)
      SMEAR_SECTION_MAP_FILE="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-interp-file)
      SMEAR_INTERP_FILE="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-data-tree)
      SMEAR_DATA_TREE="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-sim-tree)
      SMEAR_SIM_TREE="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-nx)
      SMEAR_NX="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-ny)
      SMEAR_NY="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-x-min)
      SMEAR_X_MIN="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-x-max)
      SMEAR_X_MAX="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-y-min)
      SMEAR_Y_MIN="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-y-max)
      SMEAR_Y_MAX="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-overlap)
      SMEAR_OVERLAP="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-nsmear)
      SMEAR_NSMEAR="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-sim-excl-input)
      SMEAR_SIM_EXCL_INPUT="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-sim-sidis-input)
      SMEAR_SIM_SIDIS_INPUT="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --smearing-sim-delta-input)
      SMEAR_SIM_DELTA_INPUT="$2"
      RUN_SMEARING_STAGE="yes"
      shift 2
      ;;
    --run-xsec)
      RUN_XSEC_STAGE="yes"
      shift
      ;;
    --xsec-kin)
      XSEC_KINS+=("$2")
      RUN_XSEC_STAGE="yes"
      shift 2
      ;;
    --xsec-script)
      XSEC_SCRIPT="$2"
      RUN_XSEC_STAGE="yes"
      shift 2
      ;;
    --xsec-target)
      XSEC_TARGET="$2"
      RUN_XSEC_STAGE="yes"
      shift 2
      ;;
    --xsec-root-dir)
      XSEC_ROOT_DIR="$2"
      RUN_XSEC_STAGE="yes"
      shift 2
      ;;
    --xsec-data-file)
      XSEC_DATA_FILE="$2"
      RUN_XSEC_STAGE="yes"
      shift 2
      ;;
    --xsec-sim-file)
      XSEC_SIM_FILE="$2"
      RUN_XSEC_STAGE="yes"
      shift 2
      ;;
    --xsec-out-dir)
      XSEC_OUT_DIR="$2"
      RUN_XSEC_STAGE="yes"
      shift 2
      ;;
    --xsec-out-root)
      XSEC_OUT_ROOT="$2"
      RUN_XSEC_STAGE="yes"
      shift 2
      ;;
    --xsec-out-csv)
      XSEC_OUT_CSV="$2"
      RUN_XSEC_STAGE="yes"
      shift 2
      ;;
    --xsec-out-slice-csv)
      XSEC_OUT_SLICE_CSV="$2"
      RUN_XSEC_STAGE="yes"
      shift 2
      ;;
    --xsec-all-plots-pdf)
      XSEC_ALL_PLOTS_PDF="$2"
      RUN_XSEC_STAGE="yes"
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
      shift
      if [[ $# -eq 0 || "$1" == --* ]]; then
        echo "--run requires at least one run number." >&2
        exit 1
      fi
      while [[ $# -gt 0 && "$1" != --* ]]; do
        IFS=',' read -r -a RUN_TOKENS <<< "$1"
        for run_token in "${RUN_TOKENS[@]}"; do
          run_token="$(trim_ws "${run_token}")"
          [[ -z "${run_token}" ]] && continue
          if ! [[ "${run_token}" =~ ^[0-9]+$ ]]; then
            echo "Invalid --run value: ${run_token}" >&2
            exit 1
          fi
          RUN_FILTERS+=("${run_token}")
        done
        shift
      done
      ;;
    --gevnum-cut)
      USE_GEVNUM_CUT="$2"
      shift 2
      ;;
    --jobs)
      JOBS="$2"
      shift 2
      ;;
    --timeout)
      TIMEOUT_SEC="$2"
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

if [[ ${#RUN_FILTERS[@]} -gt 0 ]]; then
  UNIQUE_RUN_FILTERS=()
  declare -A RUN_FILTER_SEEN=()
  for run in "${RUN_FILTERS[@]}"; do
    if [[ -z "${RUN_FILTER_SEEN[${run}]+x}" ]]; then
      RUN_FILTER_SEEN["${run}"]=1
      UNIQUE_RUN_FILTERS+=("${run}")
    fi
  done
  RUN_FILTERS=("${UNIQUE_RUN_FILTERS[@]}")
  RUN_FILTER="${RUN_FILTERS[*]}"
fi

if [[ "${ACCEPTANCE_CUTS_CONFIG}" != /* ]]; then
  ACCEPTANCE_CUTS_CONFIG="${REPO_ROOT}/${ACCEPTANCE_CUTS_CONFIG}"
fi

if [[ ! -f "${CONFIG_CSV}" ]]; then
  echo "Config CSV not found: ${CONFIG_CSV}" >&2
  exit 1
fi
if [[ ! -f "${ROOT_MACRO}" ]]; then
  echo "Main analysis macro not found: ${ROOT_MACRO}" >&2
  exit 1
fi
if [[ ! -f "${ACCEPTANCE_CUTS_IMPL}" ]]; then
  echo "Acceptance-cuts implementation not found: ${ACCEPTANCE_CUTS_IMPL}" >&2
  exit 1
fi
if [[ ! -f "${ACCEPTANCE_CUTS_CONFIG}" ]]; then
  echo "Acceptance-cuts config not found: ${ACCEPTANCE_CUTS_CONFIG}" >&2
  exit 1
fi
if ! command -v "${ROOT_CMD}" >/dev/null 2>&1; then
  echo "ROOT executable not found in PATH (${ROOT_CMD})." >&2
  exit 1
fi
if [[ "${COMBINE_AFTER_RUN}" == "yes" ]] && ! command -v "${PYTHON_CMD}" >/dev/null 2>&1; then
  echo "Python executable not found in PATH (${PYTHON_CMD})." >&2
  exit 1
fi
if [[ "${COMBINE_AFTER_RUN}" == "yes" && ! -f "${COMBINE_SCRIPT}" ]]; then
  echo "Combine script not found: ${COMBINE_SCRIPT}" >&2
  exit 1
fi

if [[ "${RUN_SMEARING_STAGE}" == "yes" ]]; then
  if [[ "${SMEAR_SCRIPT}" != /* ]]; then
    SMEAR_SCRIPT="${REPO_ROOT}/${SMEAR_SCRIPT}"
  fi
  if [[ ! -f "${SMEAR_SCRIPT}" ]]; then
    echo "Smearing pipeline script not found: ${SMEAR_SCRIPT}" >&2
    exit 1
  fi
  if [[ ! -x "${SMEAR_SCRIPT}" ]]; then
    echo "Smearing script is not executable; attempting to run via bash wrapper." >&2
  fi
fi

if [[ "${RUN_XSEC_STAGE}" == "yes" ]]; then
  if [[ "${XSEC_SCRIPT}" != /* ]]; then
    XSEC_SCRIPT="${REPO_ROOT}/${XSEC_SCRIPT}"
  fi
  if [[ ! -f "${XSEC_SCRIPT}" ]]; then
    echo "Xsec pipeline script not found: ${XSEC_SCRIPT}" >&2
    exit 1
  fi
  if [[ ! -x "${XSEC_SCRIPT}" ]]; then
    echo "Xsec script is not executable; attempting to run via bash wrapper." >&2
  fi
fi

if [[ ${ALL_KINS} -eq 1 && ${#KIN_LIST[@]} -gt 0 ]]; then
  echo "Use either --all-kins or --kin, not both." >&2
  exit 1
fi

SELECTED_KINS=()
if [[ ${ALL_KINS} -eq 1 ]]; then
  mapfile -t SELECTED_KINS < <(get_all_kins)
elif [[ ${#KIN_LIST[@]} -gt 0 ]]; then
  SELECTED_KINS=("${KIN_LIST[@]}")
elif [[ ${#RUN_FILTERS[@]} -gt 0 ]]; then
  declare -A SELECTED_KIN_SEEN=()
  for run in "${RUN_FILTERS[@]}"; do
    mapfile -t RUN_KINS < <(get_kins_for_run "${run}")
    if [[ ${#RUN_KINS[@]} -eq 0 ]]; then
      echo "No Kin_old entry found in config for run ${run}." >&2
      exit 1
    fi
    for kin in "${RUN_KINS[@]}"; do
      if [[ -z "${SELECTED_KIN_SEEN[${kin}]+x}" ]]; then
        SELECTED_KIN_SEEN["${kin}"]=1
        SELECTED_KINS+=("${kin}")
      fi
    done
  done
  if [[ ${#SELECTED_KINS[@]} -eq 0 ]]; then
    echo "No Kin_old entries found in config for run filter: ${RUN_FILTER}." >&2
    exit 1
  fi
  echo "[select] Auto-selected kinematic setting(s) for run(s) ${RUN_FILTER}: ${SELECTED_KINS[*]}"
else
  mapfile -t AVAILABLE_KINS < <(get_all_kins)
  if [[ ! -t 0 ]]; then
    echo "[select] No --kin/--all-kins and non-interactive shell; defaulting to all kinematics."
    SELECTED_KINS=("${AVAILABLE_KINS[@]}")
  else
    prompt_for_kinematics "${AVAILABLE_KINS[@]}"
  fi
fi

if [[ ${#SELECTED_KINS[@]} -eq 0 ]]; then
  echo "No kinematic settings selected." >&2
  exit 1
fi

SMEAR_SELECTED_KINS=()
if [[ "${RUN_SMEARING_STAGE}" == "yes" ]]; then
  if [[ -z "${SMEAR_TARGET}" ]]; then
    SMEAR_TARGET="${COMBINE_TARGET}"
  fi
  if [[ -z "${SMEAR_MODE_HINT}" ]]; then
    SMEAR_MODE_HINT="${MODE}"
  fi
  if [[ ${#SMEAR_KINS[@]} -gt 0 ]]; then
    SMEAR_SELECTED_KINS=("${SMEAR_KINS[@]}")
  else
    SMEAR_SELECTED_KINS=("${SELECTED_KINS[@]}")
  fi
fi

XSEC_SELECTED_KINS=()
if [[ "${RUN_XSEC_STAGE}" == "yes" ]]; then
  if [[ -z "${XSEC_TARGET}" ]]; then
    XSEC_TARGET="${COMBINE_TARGET}"
  fi
  if [[ ${#XSEC_KINS[@]} -gt 0 ]]; then
    XSEC_SELECTED_KINS=("${XSEC_KINS[@]}")
  else
    XSEC_SELECTED_KINS=("${SELECTED_KINS[@]}")
  fi
fi

if [[ -z "${INPUT_DIR}" ]]; then
  if [[ -n "${SOURCE}" ]]; then
    case "${SOURCE}" in
      updated)
        INPUT_DIR="${UPDATED_DIR}"
        [[ "${MODE}" == "auto" ]] && MODE="hcana"
        ;;
      production)
        INPUT_DIR="${PRODUCTION_DIR}"
        [[ "${MODE}" == "auto" ]] && MODE="hcana"
        ;;
      waveform)
        INPUT_DIR="${WAVEFORM_DIR}"
        [[ "${MODE}" == "auto" ]] && MODE="waveform"
        ;;
      *)
        echo "Invalid --source value: ${SOURCE}" >&2
        exit 1
        ;;
    esac
  else
    echo "[source] No --source/--input-dir provided; will auto-probe per run in order: waveform -> updated -> production"
  fi
fi

if [[ -n "${INPUT_DIR}" && ! -d "${INPUT_DIR}" ]]; then
  echo "Input directory does not exist: ${INPUT_DIR}" >&2
  exit 1
fi

if [[ -z "${INPUT_DIR}" && -n "${SOURCE}" ]]; then
  # Defensive check for forced source mode.
  case "${SOURCE}" in
    waveform)
      [[ -d "${WAVEFORM_DIR}" ]] || { echo "Waveform directory does not exist: ${WAVEFORM_DIR}" >&2; exit 1; }
      ;;
    updated)
      [[ -d "${UPDATED_DIR}" ]] || { echo "Updated directory does not exist: ${UPDATED_DIR}" >&2; exit 1; }
      ;;
    production)
      [[ -d "${PRODUCTION_DIR}" ]] || { echo "Production directory does not exist: ${PRODUCTION_DIR}" >&2; exit 1; }
      ;;
  esac
fi

if ! [[ "${JOBS}" =~ ^[0-9]+$ ]] || [[ "${JOBS}" -lt 1 ]]; then
  echo "Invalid --jobs value: ${JOBS}" >&2
  exit 1
fi
if ! [[ "${TIMEOUT_SEC}" =~ ^[0-9]+$ ]]; then
  echo "Invalid --timeout value: ${TIMEOUT_SEC}" >&2
  exit 1
fi

USE_GEVNUM_CUT="$(echo "${USE_GEVNUM_CUT}" | tr '[:upper:]' '[:lower:]')"
case "${USE_GEVNUM_CUT}" in
  yes|y|1|true|use|on)
    USE_GEVNUM_CUT="yes"
    ;;
  no|n|0|false|skip|off)
    USE_GEVNUM_CUT="no"
    ;;
  ask|prompt|"")
    USE_GEVNUM_CUT="ask"
    ;;
  *)
    echo "Invalid --gevnum-cut value: ${USE_GEVNUM_CUT}" >&2
    echo "Valid values: ask, yes, no" >&2
    exit 1
    ;;
esac

if [[ "${USE_GEVNUM_CUT}" == "ask" ]]; then
  if [[ -t 0 ]]; then
    prompt_for_gevnum_cut
  else
    USE_GEVNUM_CUT="yes"
    echo "[select] Non-interactive shell; defaulting g.evnum cut to 'yes'. Override with --gevnum-cut no."
  fi
fi

KIN_PATTERN="$(join_by_pipe "${SELECTED_KINS[@]}")"
TYPE_ITEMS=()
IFS=',' read -r -a TYPE_ITEMS <<< "${TYPES_CSV}"
TYPE_PATTERN="$(join_by_pipe "${TYPE_ITEMS[@]}")"
RUN_PATTERN="$(join_by_pipe "${RUN_FILTERS[@]}")"

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT
mkdir -p "${TMP_DIR}/jobs"

JOB_LIST="${TMP_DIR}/jobs.txt"

awk -F',' -v kin_pattern="${KIN_PATTERN}" -v type_pattern="${TYPE_PATTERN}" -v run_pattern="${RUN_PATTERN}" '
  function trim(s) {
    gsub(/\r/, "", s)
    sub(/^[[:space:]]+/, "", s)
    sub(/[[:space:]]+$/, "", s)
    gsub(/^"|"$/, "", s)
    return s
  }

  BEGIN {
    n_kin = split(kin_pattern, kin_parts, "|")
    for (i=1; i<=n_kin; ++i) {
      kin = trim(kin_parts[i])
      if (kin != "") wanted_kin[kin] = 1
    }

    n_type = split(type_pattern, type_parts, "|")
    for (i=1; i<=n_type; ++i) {
      typ = tolower(trim(type_parts[i]))
      if (typ != "") wanted_type[typ] = 1
    }

    n_run = split(run_pattern, run_parts, "|")
    for (i=1; i<=n_run; ++i) {
      run = trim(run_parts[i])
      if (run != "") wanted_run[run] = 1
    }
    use_run_filter = (length(run_pattern) > 0)
  }

  NR == 1 {
    for (i=1; i<=NF; ++i) {
      h = trim($i)
      if (h == "run_number") run_col = i
      if (h == "Kin_old") kin_col = i
      if (h == "Type") type_col = i
    }
    next
  }

  NR > 1 {
    if (!(kin_col > 0 && run_col > 0 && type_col > 0)) next

    kin  = trim($kin_col)
    run  = trim($run_col)
    type = tolower(trim($type_col))

    if (!(kin in wanted_kin)) next
    if (!(type in wanted_type)) next
    if (run !~ /^[0-9]+$/) next
    if (use_run_filter && !(run in wanted_run)) next

    key = kin "," run
    if (!(key in seen)) {
      print key
      seen[key] = 1
    }
  }
' "${CONFIG_CSV}" | sort > "${JOB_LIST}"

TOTAL_JOBS="$(wc -l < "${JOB_LIST}")"
if [[ "${TOTAL_JOBS}" -eq 0 ]]; then
  echo "No runs matched selected kinematics/types." >&2
  exit 1
fi

EFFECTIVE_JOBS="${JOBS}"
if [[ "${EFFECTIVE_JOBS}" -gt "${TOTAL_JOBS}" ]]; then
  EFFECTIVE_JOBS="${TOTAL_JOBS}"
fi

echo "==========================================================================="
echo "NPS π0 Unified Parallel Analysis"
echo "Start time:      $(date)"
echo "Kinematics:      ${SELECTED_KINS[*]}"
if [[ -n "${INPUT_DIR}" ]]; then
  echo "Input directory: ${INPUT_DIR}"
else
  echo "Input directory: auto-probe per run (waveform -> updated -> production)"
fi
echo "Mode:            ${MODE}"
echo "Acceptance cuts: ${ACCEPTANCE_CUTS_CONFIG}"
echo "g.evnum cut:      ${USE_GEVNUM_CUT}"
echo "Type filter:     ${TYPES_CSV}"
echo "Combine step:    ${COMBINE_AFTER_RUN} (target=${COMBINE_TARGET})"
if [[ "${RUN_SMEARING_STAGE}" == "yes" ]]; then
  echo "Smearing stage:  yes (target=${SMEAR_TARGET})"
  echo "Smearing mode:   ${SMEAR_MODE_HINT}"
  echo "Smearing script: ${SMEAR_SCRIPT}"
  echo "Smearing kins:   ${SMEAR_SELECTED_KINS[*]}"
else
  echo "Smearing stage:  no"
fi
if [[ "${RUN_XSEC_STAGE}" == "yes" ]]; then
  echo "Xsec stage:      yes (target=${XSEC_TARGET})"
  echo "Xsec script:     ${XSEC_SCRIPT}"
  echo "Xsec kins:       ${XSEC_SELECTED_KINS[*]}"
else
  echo "Xsec stage:      no"
fi
echo "Run filter:      ${RUN_FILTER:-<none>}"
echo "Parallel jobs:   ${EFFECTIVE_JOBS}/${TOTAL_JOBS}"
echo "Output base:     ${OUTPUT_BASE}"
echo "==========================================================================="

PROGRESS_FILE="${TMP_DIR}/progress.done"
: > "${PROGRESS_FILE}"

export ROOT_CMD ROOT_MACRO ACCEPTANCE_CUTS_IMPL CONFIG_CSV OUTPUT_BASE INPUT_DIR MODE TYPES_CSV TIMEOUT_SEC TMP_DIR PROGRESS_FILE CSV_HEADER USE_GEVNUM_CUT
export UPDATED_DIR PRODUCTION_DIR WAVEFORM_DIR SOURCE
export ACCEPTANCE_CUTS_CONFIG

xargs -P "${EFFECTIVE_JOBS}" -I {} bash -c '
  line="$1"
  IFS="," read -r kin run <<< "${line}"
  safe_kin="$(echo "${kin}" | sed "s/[^[:alnum:]_-]/_/g")"

  job_dir="${TMP_DIR}/jobs/${safe_kin}"
  mkdir -p "${job_dir}"
  run_log_dir="${OUTPUT_BASE}/${safe_kin}/logs"
  mkdir -p "${run_log_dir}"

  run_log="${run_log_dir}/analysis_main_run${run}.log"
  status_file="${job_dir}/status_${run}.txt"
  row_csv_copy="${job_dir}/run_${run}.csv"
  resolved_input_dir="${INPUT_DIR}"
  resolved_mode="${MODE}"

  # Resolve input ROOT location for this run.
  # Priority: explicit --input-dir, forced --source, else auto-probe
  # order waveform -> updated -> production.
  if [[ -z "${resolved_input_dir}" ]]; then
    case "${SOURCE}" in
      waveform)
        resolved_input_dir="${WAVEFORM_DIR}"
        [[ "${resolved_mode}" == "auto" ]] && resolved_mode="waveform"
        ;;
      updated)
        resolved_input_dir="${UPDATED_DIR}"
        [[ "${resolved_mode}" == "auto" ]] && resolved_mode="hcana"
        ;;
      production)
        resolved_input_dir="${PRODUCTION_DIR}"
        [[ "${resolved_mode}" == "auto" ]] && resolved_mode="hcana"
        ;;
      "")
        if compgen -G "${WAVEFORM_DIR%/}/nps_production_${run}_*_wf_calib.root" >/dev/null; then
          resolved_input_dir="${WAVEFORM_DIR}"
          [[ "${resolved_mode}" == "auto" ]] && resolved_mode="waveform"
        elif compgen -G "${UPDATED_DIR%/}/skim_run${run}.root" >/dev/null || \
             compgen -G "${UPDATED_DIR%/}/nps_hms_coin_${run}_*_1_-1.root" >/dev/null; then
          resolved_input_dir="${UPDATED_DIR}"
          [[ "${resolved_mode}" == "auto" ]] && resolved_mode="hcana"
        elif compgen -G "${PRODUCTION_DIR%/}/skim_run${run}.root" >/dev/null || \
             compgen -G "${PRODUCTION_DIR%/}/nps_hms_coin_${run}_*_1_-1.root" >/dev/null; then
          resolved_input_dir="${PRODUCTION_DIR}"
          [[ "${resolved_mode}" == "auto" ]] && resolved_mode="hcana"
        else
          {
            echo "[ERROR] No input ROOT files found for run ${run} in probe order:"
            echo "        waveform:   ${WAVEFORM_DIR}"
            echo "        updated:    ${UPDATED_DIR}"
            echo "        production: ${PRODUCTION_DIR}"
          } > "${run_log}" 2>&1
          echo "${kin},${run},90,${row_csv_copy}" > "${status_file}"
          echo "${kin},${run},90" >> "${PROGRESS_FILE}"
          exit 90
        fi
        ;;
      *)
        {
          echo "[ERROR] Invalid SOURCE value: ${SOURCE}"
        } > "${run_log}" 2>&1
        echo "${kin},${run},91,${row_csv_copy}" > "${status_file}"
        echo "${kin},${run},91" >> "${PROGRESS_FILE}"
        exit 91
        ;;
    esac
  fi

  export NPS_KIN="${kin}"
  export NPS_MODE="${resolved_mode}"
  export NPS_CONFIG_CSV="${CONFIG_CSV}"
  export NPS_ACCEPTANCE_CUTS_CONFIG="${ACCEPTANCE_CUTS_CONFIG}"
  export NPS_OUTPUT_BASE="${OUTPUT_BASE}"
  export NPS_INPUT_DIR="${resolved_input_dir}"
  export NPS_TYPES="${TYPES_CSV}"
  export NPS_RUN="${run}"
  export NPS_USE_GEVNUM_CUT="${USE_GEVNUM_CUT}"

  status=0
  {
    echo "[INFO] run=${run} kin=${kin} input_dir=${resolved_input_dir} mode=${resolved_mode} gevnum_cut=${USE_GEVNUM_CUT}"
    if [[ "${TIMEOUT_SEC}" -gt 0 ]]; then
      timeout "${TIMEOUT_SEC}" "${ROOT_CMD}" -l -b -e ".L ${ACCEPTANCE_CUTS_IMPL}+" -q "${ROOT_MACRO}()"
    else
      "${ROOT_CMD}" -l -b -e ".L ${ACCEPTANCE_CUTS_IMPL}+" -q "${ROOT_MACRO}()"
    fi
  } > "${run_log}" 2>&1 || status=$?

  csv_src=""
  if grep -q "\[CSV_WRITTEN\]" "${run_log}" 2>/dev/null; then
    csv_src="$(grep "\[CSV_WRITTEN\]" "${run_log}" | tail -1 | awk "{print \$NF}")"
    if [[ -f "${csv_src}" ]]; then
      cp "${csv_src}" "${row_csv_copy}"
    fi
  fi

  if [[ "${status}" -eq 0 && ! -s "${row_csv_copy}" ]]; then
    {
      echo "[ERROR] Analysis exited with status 0 but produced no per-run CSV row."
      echo "[ERROR] Expected copied row file: ${row_csv_copy}"
      echo "[ERROR] CSV source marker: ${csv_src:-<none>}"
    } >> "${run_log}" 2>&1
    status=92
  fi

  echo "${kin},${run},${status},${row_csv_copy}" > "${status_file}"
  echo "${kin},${run},${status}" >> "${PROGRESS_FILE}"
  exit "${status}"
' _ {} < "${JOB_LIST}" || true

FAILED=0

declare -A SMEAR_KIN_SET=()
if [[ "${RUN_SMEARING_STAGE}" == "yes" ]]; then
  for kin in "${SMEAR_SELECTED_KINS[@]}"; do
    SMEAR_KIN_SET["${kin}"]=1
  done
fi

declare -A XSEC_KIN_SET=()
if [[ "${RUN_XSEC_STAGE}" == "yes" ]]; then
  for kin in "${XSEC_SELECTED_KINS[@]}"; do
    XSEC_KIN_SET["${kin}"]=1
  done
fi

while IFS= read -r status_file; do
  [[ -f "${status_file}" ]] || continue
  IFS=',' read -r kin run status row_csv < "${status_file}"
  if [[ "${status}" -ne 0 ]]; then
    echo "[fail] kin=${kin} run=${run} status=${status}"
    FAILED=$((FAILED + 1))
  fi
done < <(find "${TMP_DIR}/jobs" -type f -name 'status_*.txt' | sort)

for kin in "${SELECTED_KINS[@]}"; do
  safe_kin="$(sanitize_name "${kin}")"
  summary_dir="${OUTPUT_BASE}/${safe_kin}/summary"
  mkdir -p "${summary_dir}"
  summary_csv="${summary_dir}/summary_all_runs.csv"

  echo "${CSV_HEADER}" > "${summary_csv}"
  mapfile -t rows < <(find "${TMP_DIR}/jobs/${safe_kin}" -type f -name 'run_*.csv' | sort)
  for row_file in "${rows[@]}"; do
    tail -n 1 "${row_file}" >> "${summary_csv}"
  done

  if [[ ${#rows[@]} -eq 0 ]]; then
    echo "[fail] No per-run CSV rows found for kin=${kin}" >&2
    FAILED=$((FAILED + 1))
    continue
  else
    echo "[merge] Rebuilt ${summary_csv} from ${#rows[@]} run row files"
  fi

  combine_cut_debug_pdfs "${kin}"

  if [[ "${COMBINE_AFTER_RUN}" == "yes" ]]; then
    kin_root_dir="${OUTPUT_BASE}/${safe_kin}/root"
    combine_log="${OUTPUT_BASE}/${safe_kin}/logs/combine_${safe_kin}.log"
    echo "[combine] Running combine stage for kin=${kin} (target=${COMBINE_TARGET}, root_dir=${kin_root_dir})"
    status=0
    "${PYTHON_CMD}" "${COMBINE_SCRIPT}" \
      --kin "${kin}" \
      --config "${CONFIG_CSV}" \
      --output-base "${OUTPUT_BASE}" \
      --root-dir "${kin_root_dir}" \
      --target "${COMBINE_TARGET}" \
      --types "${TYPES_CSV}" \
      # --fp-debug-plots \
      > "${combine_log}" 2>&1 || status=$?

    if [[ "${status}" -ne 0 ]]; then
      echo "[fail] combine step failed for kin=${kin} status=${status}; see ${combine_log}" >&2
      FAILED=$((FAILED + 1))
    else
      echo "[combine] Completed kin=${kin}; log=${combine_log}"
    fi
  fi

  if [[ "${RUN_SMEARING_STAGE}" == "yes" && -n "${SMEAR_KIN_SET[${kin}]+x}" ]]; then
    smear_log="${OUTPUT_BASE}/${safe_kin}/logs/smearing_${safe_kin}.log"
    echo "[smear] Running simulation smearing for kin=${kin} (target=${SMEAR_TARGET})"

    smear_cmd=("${SMEAR_SCRIPT}"
      --kin "${kin}"
      --target "${SMEAR_TARGET}"
      --output-base "${OUTPUT_BASE}"
      --mode "${SMEAR_MODE_HINT}"
      --acceptance-config "${ACCEPTANCE_CUTS_CONFIG}")

    [[ -n "${SMEAR_ROOT_DIR}" ]] && smear_cmd+=(--root-dir "${SMEAR_ROOT_DIR}")
    [[ -n "${SMEAR_COMBINED_FILE}" ]] && smear_cmd+=(--combined-file "${SMEAR_COMBINED_FILE}")
    [[ -n "${SMEAR_SIM_FILE}" ]] && smear_cmd+=(--sim-file "${SMEAR_SIM_FILE}")
    [[ -n "${SMEAR_SIM_SMEARED_FILE}" ]] && smear_cmd+=(--sim-smeared-file "${SMEAR_SIM_SMEARED_FILE}")
    [[ -n "${SMEAR_OUT_FILE}" ]] && smear_cmd+=(--out-file "${SMEAR_OUT_FILE}")
    [[ -n "${SMEAR_SECTION_MAP_FILE}" ]] && smear_cmd+=(--section-map-file "${SMEAR_SECTION_MAP_FILE}")
    [[ -n "${SMEAR_INTERP_FILE}" ]] && smear_cmd+=(--interp-file "${SMEAR_INTERP_FILE}")
    [[ -n "${SMEAR_DATA_TREE}" ]] && smear_cmd+=(--data-tree "${SMEAR_DATA_TREE}")
    [[ -n "${SMEAR_SIM_TREE}" ]] && smear_cmd+=(--sim-tree "${SMEAR_SIM_TREE}")
    [[ -n "${SMEAR_NX}" ]] && smear_cmd+=(--nx "${SMEAR_NX}")
    [[ -n "${SMEAR_NY}" ]] && smear_cmd+=(--ny "${SMEAR_NY}")
    [[ -n "${SMEAR_X_MIN}" ]] && smear_cmd+=(--x-min "${SMEAR_X_MIN}")
    [[ -n "${SMEAR_X_MAX}" ]] && smear_cmd+=(--x-max "${SMEAR_X_MAX}")
    [[ -n "${SMEAR_Y_MIN}" ]] && smear_cmd+=(--y-min "${SMEAR_Y_MIN}")
    [[ -n "${SMEAR_Y_MAX}" ]] && smear_cmd+=(--y-max "${SMEAR_Y_MAX}")
    [[ -n "${SMEAR_OVERLAP}" ]] && smear_cmd+=(--overlap "${SMEAR_OVERLAP}")
    [[ -n "${SMEAR_NSMEAR}" ]] && smear_cmd+=(--nsmear "${SMEAR_NSMEAR}")
    [[ -n "${SMEAR_SIM_EXCL_INPUT}" ]] && smear_cmd+=(--sim-excl-input "${SMEAR_SIM_EXCL_INPUT}")
    [[ -n "${SMEAR_SIM_SIDIS_INPUT}" ]] && smear_cmd+=(--sim-sidis-input "${SMEAR_SIM_SIDIS_INPUT}")
    [[ -n "${SMEAR_SIM_DELTA_INPUT}" ]] && smear_cmd+=(--sim-delta-input "${SMEAR_SIM_DELTA_INPUT}")

    status=0
    bash "${smear_cmd[@]}" > "${smear_log}" 2>&1 || status=$?
    if [[ "${status}" -ne 0 ]]; then
      echo "[fail] smearing stage failed for kin=${kin} status=${status}; see ${smear_log}" >&2
      FAILED=$((FAILED + 1))
    else
      echo "[smear] Completed kin=${kin}; log=${smear_log}"
    fi
  fi

  if [[ "${RUN_XSEC_STAGE}" == "yes" && -n "${XSEC_KIN_SET[${kin}]+x}" ]]; then
    xsec_log="${OUTPUT_BASE}/${safe_kin}/logs/xsec_${safe_kin}.log"
    echo "[xsec] Running xsec extraction for kin=${kin} (target=${XSEC_TARGET})"

    xsec_cmd=("${XSEC_SCRIPT}"
      --kin "${kin}"
      --target "${XSEC_TARGET}"
      --output-base "${OUTPUT_BASE}")

    [[ -n "${XSEC_ROOT_DIR}" ]] && xsec_cmd+=(--root-dir "${XSEC_ROOT_DIR}")
    [[ -n "${XSEC_DATA_FILE}" ]] && xsec_cmd+=(--data-file "${XSEC_DATA_FILE}")
    [[ -n "${XSEC_SIM_FILE}" ]] && xsec_cmd+=(--sim-file "${XSEC_SIM_FILE}")
    [[ -n "${XSEC_OUT_DIR}" ]] && xsec_cmd+=(--out-dir "${XSEC_OUT_DIR}")
    [[ -n "${XSEC_OUT_ROOT}" ]] && xsec_cmd+=(--out-root "${XSEC_OUT_ROOT}")
    [[ -n "${XSEC_OUT_CSV}" ]] && xsec_cmd+=(--out-csv "${XSEC_OUT_CSV}")
    [[ -n "${XSEC_OUT_SLICE_CSV}" ]] && xsec_cmd+=(--out-slice-csv "${XSEC_OUT_SLICE_CSV}")
    [[ -n "${XSEC_ALL_PLOTS_PDF}" ]] && xsec_cmd+=(--all-plots-pdf "${XSEC_ALL_PLOTS_PDF}")

    status=0
    bash "${xsec_cmd[@]}" > "${xsec_log}" 2>&1 || status=$?
    if [[ "${status}" -ne 0 ]]; then
      echo "[fail] xsec stage failed for kin=${kin} status=${status}; see ${xsec_log}" >&2
      FAILED=$((FAILED + 1))
    else
      echo "[xsec] Completed kin=${kin}; log=${xsec_log}"
    fi
  fi
done

echo "==========================================================================="
echo "Unified parallel analysis complete at $(date)"
echo "Total jobs: ${TOTAL_JOBS}, failed: ${FAILED}"
echo "==========================================================================="

if [[ "${FAILED}" -ne 0 ]]; then
  echo "Some runs failed. Check per-run logs under ${OUTPUT_BASE}/<kin>/logs" >&2
  exit 1
fi

echo "All runs completed successfully."
