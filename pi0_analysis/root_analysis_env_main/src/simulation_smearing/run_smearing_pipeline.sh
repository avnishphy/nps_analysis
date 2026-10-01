#!/usr/bin/env bash
# ============================================================================
# End-to-end NPS smearing pipeline
# 1) Resolve and validate the raw SIMC/Geant production inputs
# 2) Generate a fresh nominal simulation into current-run staging
# 3) Fit section coefficients and generate interpolated maps from that exact file
# 4) Replay the same raw production in smeared-only mode
# 5) Validate everything, then publish and archive the completed run atomically per file
# ============================================================================

# ./src/simulation_smearing/run_smearing_pipeline.sh --kin KinC_x60_4b --target LH2 --combined-file output/KinC_x60_4b_pepsi_15gev/KinC_x60_4b/root/combined_branches_LH2.root --output-base output/KinC_x60_4b_pepsi_15gev/ --nx 7 --ny 8 --overlap 0.1

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

ROOT_CMD="${ROOT_CMD:-root}"
CXX_CMD="${CXX:-g++}"
RUN_TAG="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)_$$}"
export RUN_TAG

SMEAR_SRC="${SCRIPT_DIR}/nps_sim_smearing_new.C"
SIMC_MACRO="${SCRIPT_DIR}/simc_pi0_analysis.C"
PUBLISH_HELPER="${SCRIPT_DIR}/publish_smearing_run.py"

KIN="${NPS_SMEAR_KIN:-x60_4b}"
TARGET="${NPS_SMEAR_TARGET:-LH2}"
OUTPUT_BASE="${NPS_OUTPUT_BASE:-${REPO_ROOT}/output}"

ROOT_DIR=""
COMBINED_FILE=""
SIM_FILE=""
SIM_SMEARED_FILE=""
OUT_FILE=""
SECTION_MAP_FILE=""
INTERP_FILE=""
INPUT_PREVIEW_PDF="${NPS_SMEAR_INPUT_PREVIEW_PDF:-}"

DATA_TREE="${NPS_SMEAR_DATA_TREE:-physics}"
SIM_TREE="${NPS_SMEAR_SIM_TREE:-simulation}"

NX="${NPS_SMEAR_NX:-11}"
NY="${NPS_SMEAR_NY:-16}"
X_MIN="${NPS_SMEAR_X_MIN:--24}"
X_MAX="${NPS_SMEAR_X_MAX:-28}"
Y_MIN="${NPS_SMEAR_Y_MIN:--34}"
Y_MAX="${NPS_SMEAR_Y_MAX:-34}"
OVERLAP="${NPS_SMEAR_OVERLAP:-0.0}"
NSMEAR="${NPS_SMEAR_NSMEAR:-80}"
SMEAR_RANDOM_SEED="${NPS_SMEAR_RANDOM_SEED:-42}"
NPS_Z_NPS_CM_ARG="${NPS_Z_NPS_CM:-${NPS_SMEAR_Z_NPS_CM:-}}"
NPS_THETA_NPS_DEG_ARG="${NPS_THETA_NPS_DEG:-${NPS_SMEAR_THETA_NPS_DEG:-}}"
NPS_ANGLE_DEG_ARG="${NPS_NPS_THETA_DEG:-${NPS_SMEAR_NPS_ANGLE_DEG:-}}"
BEAM_ENERGY_GEV_ARG="${NPS_EBEAM:-${NPS_SMEAR_BEAM_ENERGY_GEV:-}}"

SIM_INPUT_EXCLUSIVE="${NPS_SIMC_EXCLUSIVE_INPUT:-${NPS_SIMC_INPUT_FILE_EXCLUSIVE:-}}"
SIM_INPUT_SIDIS="${NPS_SIMC_SIDIS_INPUT:-${NPS_SIMC_INPUT_FILE_SIDIS:-}}"
SIM_INPUT_DELTA="${NPS_SIMC_DELTA_INPUT:-${NPS_SIMC_INPUT_FILE_DELTA:-}}"
SIMC_PRODUCTION_DIR="${NPS_SIMC_PRODUCTION_DIR:-${REPO_ROOT}/output/simc/nps_simc_20260824_135058}"
GEANT_PRODUCTION_DIR="${NPS_GEANT_PRODUCTION_DIR:-/lustre24/expphy/volatile/hallc/nps/singhav/geant4_simc/nps_geant4_20260824_221705}"
SIM_HIST_EXCLUSIVE="${NPS_SIMC_EXCLUSIVE_HIST:-}"
SIM_HIST_SIDIS="${NPS_SIMC_SIDIS_HIST:-}"
SIM_HIST_DELTA="${NPS_SIMC_DELTA_HIST:-}"
DEAD_BLOCK_CONFIG="${NPS_DEAD_BLOCK_CONFIG:-${REPO_ROOT}/config/dead_block_per_runs.csv}"
DEAD_BLOCK_RUN="${NPS_SIMC_DEAD_BLOCK_RUN:-}"
KINEMATICS_CONFIG="${NPS_SIMULATION_KINEMATICS_CONFIG:-${REPO_ROOT}/config/nps_simulation_kinematics.csv}"
MODE_HINT="${NPS_SMEAR_MODE_HINT:-${NPS_MODE:-auto}}"
ACCEPTANCE_CONFIG="${NPS_ACCEPTANCE_CUTS_CONFIG:-${REPO_ROOT}/config/acceptance_cuts.conf}"

trim_ws() {
  local s="$1"
  s="${s#${s%%[![:space:]]*}}"
  s="${s%${s##*[![:space:]]}}"
  echo "$s"
}

sanitize_name() {
  echo "$1" | sed 's/[^[:alnum:]_-]/_/g'
}

to_abs_path() {
  local p="$1"
  if [[ -z "${p}" ]]; then
    echo ""
  elif [[ "${p}" = /* ]]; then
    echo "${p}"
  else
    echo "${REPO_ROOT}/${p}"
  fi
}

derive_interp_path() {
  local p="$1"
  if [[ "${p}" == *.root ]]; then
    echo "${p%.root}_interpolated.root"
  else
    echo "${p}_interpolated.root"
  fi
}

print_help() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --kin <Kin_old>             Kinematic setting (default: x60_4b)
  --target <name>             Combined target token (default: LH2)
  --output-base <path>        Output base directory (default: repo/output)
  --root-dir <path>           Directory containing combined/sim/root smearing files
  --combined-file <path>      Combined data ROOT file (default: <root-dir>/combined_branches_<target>.root)
  --sim-file <path>           Unsmeared simulation ROOT file (default: <root-dir>/simc_pi0_analysis_output.root)
  --sim-smeared-file <path>   Smeared simulation ROOT file (default: <root-dir>/simc_pi0_analysis_output_smeared.root)
  --out-file <path>           Smearing fit output ROOT (default: <root-dir>/out_smear.root)
  --section-map-file <path>   Section CSV output path (default: dirname(<out-file>)/section_map.csv)
  --interp-file <path>        Interpolated map ROOT path (default: <out-file> with _interpolated)
  --input-preview-pdf <path>  Pre-fit data/simulation verification PDF
  --data-tree <name>          Data tree name for fitter (default: physics)
  --sim-tree <name>           Simulation tree name for fitter (default: simulation)
  --nx <int>                  Section count in x (default: 11)
  --ny <int>                  Section count in y (default: 16)
  --x-min <float>             Calorimeter x minimum (default: -24)
  --x-max <float>             Calorimeter x maximum (default: 28)
  --y-min <float>             Calorimeter y minimum (default: -34)
  --y-max <float>             Calorimeter y maximum (default: 34)
  --overlap <float>           Section overlap fraction (default: 0.0)
  --nsmear <int>              Smearing iterations/event in fit (default: 80)
  --smear-seed <int>          Record requested producer RNG seed (legacy producer RNG is unchanged)
  --beam-energy-gev <float>   Override beam energy (default: kinematics config)
  --z-nps-cm <float>          Override NPS target-to-calorimeter distance
  --theta-nps-deg <float>     Override signed NPS->Hall rotation angle
  --nps-angle-deg <float>     Override physical NPS angle; fitter uses its negative
  --mode <auto|hcana|waveform>
                               Mode hint for timing defaults in SIM producer (default: auto)
  --acceptance-config <path>  Acceptance-cuts config file (default: repo/config/acceptance_cuts.conf)
  --sim-excl-input <path>     Exclusive SIMC input (forwarded to simc_pi0_analysis.C)
  --sim-sidis-input <path>    SIDIS SIMC input (forwarded to simc_pi0_analysis.C)
  --sim-delta-input <path>    Delta pi0 SIMC input (forwarded to simc_pi0_analysis.C)
  --simc-production-dir <path>
                               SIMC production base containing channel hist files
  --geant-production-dir <path>
                               Geant production base containing excl/sidis/delta inputs
  --sim-excl-hist <path>      Exclusive SIMC hist override
  --sim-sidis-hist <path>     SIDIS SIMC hist override
  --sim-delta-hist <path>     Delta SIMC hist override
  --dead-block-config <path>  Per-run dead-block CSV
  --dead-block-run <run>      Run whose mask is encoded in output branches
  --kinematics-config <path>  CSV providing beam energy, NPS z, and theta by Kin_old
  --help                      Show this help

Environment overrides:
  ROOT_CMD, CXX,
  NPS_SMEAR_KIN, NPS_SMEAR_TARGET, NPS_OUTPUT_BASE,
  NPS_MODE, NPS_SMEAR_MODE_HINT,
  NPS_ACCEPTANCE_CUTS_CONFIG,
  NPS_SMEAR_DATA_TREE, NPS_SMEAR_SIM_TREE,
  NPS_SMEAR_NX, NPS_SMEAR_NY, NPS_SMEAR_X_MIN, NPS_SMEAR_X_MAX,
  NPS_SMEAR_Y_MIN, NPS_SMEAR_Y_MAX, NPS_SMEAR_OVERLAP, NPS_SMEAR_NSMEAR,
  NPS_SMEAR_RANDOM_SEED,
  NPS_SMEAR_INPUT_PREVIEW_PDF,
  NPS_EBEAM, NPS_SMEAR_BEAM_ENERGY_GEV,
  NPS_Z_NPS_CM, NPS_THETA_NPS_DEG, NPS_NPS_THETA_DEG,
  NPS_SMEAR_Z_NPS_CM, NPS_SMEAR_THETA_NPS_DEG, NPS_SMEAR_NPS_ANGLE_DEG,
  NPS_SIMC_EXCLUSIVE_INPUT, NPS_SIMC_SIDIS_INPUT, NPS_SIMC_DELTA_INPUT,
  NPS_SIMC_PRODUCTION_DIR, NPS_GEANT_PRODUCTION_DIR,
  NPS_SIMC_EXCLUSIVE_HIST, NPS_SIMC_SIDIS_HIST, NPS_SIMC_DELTA_HIST,
  NPS_SIMC_EXCLUSIVE_NORMFAC, NPS_SIMC_SIDIS_NORMFAC, NPS_SIMC_DELTA_NORMFAC,
  NPS_SIMC_EXCLUSIVE_NGEN, NPS_SIMC_SIDIS_NGEN, NPS_SIMC_DELTA_NGEN,
  NPS_DEAD_BLOCK_CONFIG, NPS_SIMC_DEAD_BLOCK_RUN,
  NPS_SIMULATION_KINEMATICS_CONFIG,
  legacy NPS_SIMC_INPUT_FILE_EXCLUSIVE, NPS_SIMC_INPUT_FILE_SIDIS,
  NPS_SIMC_INPUT_FILE_DELTA.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --kin)
      KIN="$2"
      shift 2
      ;;
    --target)
      TARGET="$2"
      shift 2
      ;;
    --output-base)
      OUTPUT_BASE="$2"
      shift 2
      ;;
    --root-dir)
      ROOT_DIR="$2"
      shift 2
      ;;
    --combined-file)
      COMBINED_FILE="$2"
      shift 2
      ;;
    --sim-file)
      SIM_FILE="$2"
      shift 2
      ;;
    --sim-smeared-file)
      SIM_SMEARED_FILE="$2"
      shift 2
      ;;
    --out-file)
      OUT_FILE="$2"
      shift 2
      ;;
    --section-map-file)
      SECTION_MAP_FILE="$2"
      shift 2
      ;;
    --interp-file)
      INTERP_FILE="$2"
      shift 2
      ;;
    --input-preview-pdf)
      INPUT_PREVIEW_PDF="$2"
      shift 2
      ;;
    --data-tree)
      DATA_TREE="$2"
      shift 2
      ;;
    --sim-tree)
      SIM_TREE="$2"
      shift 2
      ;;
    --nx)
      NX="$2"
      shift 2
      ;;
    --ny)
      NY="$2"
      shift 2
      ;;
    --x-min)
      X_MIN="$2"
      shift 2
      ;;
    --x-max)
      X_MAX="$2"
      shift 2
      ;;
    --y-min)
      Y_MIN="$2"
      shift 2
      ;;
    --y-max)
      Y_MAX="$2"
      shift 2
      ;;
    --overlap)
      OVERLAP="$2"
      shift 2
      ;;
    --nsmear)
      NSMEAR="$2"
      shift 2
      ;;
    --smear-seed)
      SMEAR_RANDOM_SEED="$2"
      shift 2
      ;;
    --beam-energy-gev|--beam-energy)
      BEAM_ENERGY_GEV_ARG="$2"
      shift 2
      ;;
    --z-nps-cm)
      NPS_Z_NPS_CM_ARG="$2"
      shift 2
      ;;
    --theta-nps-deg)
      NPS_THETA_NPS_DEG_ARG="$2"
      shift 2
      ;;
    --nps-angle-deg)
      NPS_ANGLE_DEG_ARG="$2"
      shift 2
      ;;
    --mode)
      MODE_HINT="$2"
      shift 2
      ;;
    --acceptance-config)
      ACCEPTANCE_CONFIG="$2"
      shift 2
      ;;
    --sim-excl-input)
      SIM_INPUT_EXCLUSIVE="$2"
      shift 2
      ;;
    --sim-sidis-input)
      SIM_INPUT_SIDIS="$2"
      shift 2
      ;;
    --sim-delta-input)
      SIM_INPUT_DELTA="$2"
      shift 2
      ;;
    --simc-production-dir)
      SIMC_PRODUCTION_DIR="$2"
      shift 2
      ;;
    --geant-production-dir)
      GEANT_PRODUCTION_DIR="$2"
      shift 2
      ;;
    --sim-excl-hist)
      SIM_HIST_EXCLUSIVE="$2"
      shift 2
      ;;
    --sim-sidis-hist)
      SIM_HIST_SIDIS="$2"
      shift 2
      ;;
    --sim-delta-hist)
      SIM_HIST_DELTA="$2"
      shift 2
      ;;
    --dead-block-config)
      DEAD_BLOCK_CONFIG="$2"
      shift 2
      ;;
    --dead-block-run)
      DEAD_BLOCK_RUN="$2"
      shift 2
      ;;
    --kinematics-config)
      KINEMATICS_CONFIG="$2"
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

KIN="$(trim_ws "${KIN}")"
TARGET="$(trim_ws "${TARGET}")"

if [[ -z "${KIN}" ]]; then
  echo "[ERROR] --kin cannot be empty" >&2
  exit 1
fi
if [[ -z "${TARGET}" ]]; then
  echo "[ERROR] --target cannot be empty" >&2
  exit 1
fi

KIN_SAFE="$(sanitize_name "${KIN}")"
TARGET_SAFE="$(sanitize_name "${TARGET}")"
OUTPUT_BASE="$(to_abs_path "${OUTPUT_BASE}")"
MODE_HINT="$(trim_ws "${MODE_HINT}")"
ACCEPTANCE_CONFIG="$(to_abs_path "${ACCEPTANCE_CONFIG}")"
SIMC_PRODUCTION_DIR="$(to_abs_path "${SIMC_PRODUCTION_DIR}")"
GEANT_PRODUCTION_DIR="$(to_abs_path "${GEANT_PRODUCTION_DIR}")"
DEAD_BLOCK_CONFIG="$(to_abs_path "${DEAD_BLOCK_CONFIG}")"
KINEMATICS_CONFIG="$(to_abs_path "${KINEMATICS_CONFIG}")"

if [[ -z "${MODE_HINT}" ]]; then
  MODE_HINT="auto"
fi

if [[ -z "${ROOT_DIR}" ]]; then
  ROOT_DIR="${OUTPUT_BASE}/${KIN_SAFE}/root"
fi

ROOT_DIR="$(to_abs_path "${ROOT_DIR}")"
if [[ -z "${COMBINED_FILE}" ]]; then
  COMBINED_FILE="${ROOT_DIR}/combined_branches_${TARGET_SAFE}.root"
fi
if [[ -z "${SIM_FILE}" ]]; then
  SIM_FILE="${ROOT_DIR}/simc_pi0_analysis_output.root"
fi
if [[ -z "${SIM_SMEARED_FILE}" ]]; then
  SIM_SMEARED_FILE="${ROOT_DIR}/simc_pi0_analysis_output_smeared.root"
fi
if [[ -z "${OUT_FILE}" ]]; then
  OUT_FILE="${ROOT_DIR}/out_smear.root"
fi

COMBINED_FILE="$(to_abs_path "${COMBINED_FILE}")"
SIM_FILE="$(to_abs_path "${SIM_FILE}")"
SIM_SMEARED_FILE="$(to_abs_path "${SIM_SMEARED_FILE}")"
OUT_FILE="$(to_abs_path "${OUT_FILE}")"

GENERATED_SECTION_MAP_FILE="$(dirname "${OUT_FILE}")/section_map.csv"
GENERATED_INTERP_FILE="$(derive_interp_path "${OUT_FILE}")"

if [[ -z "${SECTION_MAP_FILE}" ]]; then
  SECTION_MAP_FILE="${GENERATED_SECTION_MAP_FILE}"
fi
if [[ -z "${INTERP_FILE}" ]]; then
  INTERP_FILE="${GENERATED_INTERP_FILE}"
fi
if [[ -z "${INPUT_PREVIEW_PDF}" ]]; then
  INPUT_PREVIEW_PDF="${ROOT_DIR}/smearing_input_histograms_prefit.pdf"
fi

SECTION_MAP_FILE="$(to_abs_path "${SECTION_MAP_FILE}")"
INTERP_FILE="$(to_abs_path "${INTERP_FILE}")"
INPUT_PREVIEW_PDF="$(to_abs_path "${INPUT_PREVIEW_PDF}")"

if [[ -n "${SIM_INPUT_EXCLUSIVE}" ]]; then
  SIM_INPUT_EXCLUSIVE="$(to_abs_path "${SIM_INPUT_EXCLUSIVE}")"
fi
if [[ -n "${SIM_INPUT_SIDIS}" ]]; then
  SIM_INPUT_SIDIS="$(to_abs_path "${SIM_INPUT_SIDIS}")"
fi
if [[ -n "${SIM_INPUT_DELTA}" ]]; then
  SIM_INPUT_DELTA="$(to_abs_path "${SIM_INPUT_DELTA}")"
fi
if [[ -n "${SIM_HIST_EXCLUSIVE}" ]]; then
  SIM_HIST_EXCLUSIVE="$(to_abs_path "${SIM_HIST_EXCLUSIVE}")"
fi
if [[ -n "${SIM_HIST_SIDIS}" ]]; then
  SIM_HIST_SIDIS="$(to_abs_path "${SIM_HIST_SIDIS}")"
fi
if [[ -n "${SIM_HIST_DELTA}" ]]; then
  SIM_HIST_DELTA="$(to_abs_path "${SIM_HIST_DELTA}")"
fi

KIN_TOKEN="${KIN#KinC_}"
if [[ -z "${KIN_TOKEN}" || ! "${KIN_TOKEN}" =~ ^[[:alnum:]_]+$ ]]; then
  echo "[ERROR] Invalid kinematic token: ${KIN}" >&2
  exit 1
fi

SIM_HIST_DIR="${SIMC_PRODUCTION_DIR}/outfiles/simc_gfortran_updated/outfiles"
SIM_INPUT_EXCLUSIVE="${SIM_INPUT_EXCLUSIVE:-${GEANT_PRODUCTION_DIR}/excl/nps_excl_pi0_${KIN_TOKEN}_geant4.root}"
SIM_INPUT_SIDIS="${SIM_INPUT_SIDIS:-${GEANT_PRODUCTION_DIR}/sidis/nps_sidis_pi0_${KIN_TOKEN}_geant4.root}"
SIM_INPUT_DELTA="${SIM_INPUT_DELTA:-${GEANT_PRODUCTION_DIR}/delta/nps_delta_pi0_${KIN_TOKEN}_geant4.root}"
SIM_HIST_EXCLUSIVE="${SIM_HIST_EXCLUSIVE:-${SIM_HIST_DIR}/nps_excl_pi0_${KIN_TOKEN}.hist}"
SIM_HIST_SIDIS="${SIM_HIST_SIDIS:-${SIM_HIST_DIR}/nps_sidis_pi0_${KIN_TOKEN}.hist}"
SIM_HIST_DELTA="${SIM_HIST_DELTA:-${SIM_HIST_DIR}/nps_delta_pi0_${KIN_TOKEN}.hist}"

if [[ -z "${BEAM_ENERGY_GEV_ARG}" || -z "${NPS_Z_NPS_CM_ARG}" ||
      ( -z "${NPS_THETA_NPS_DEG_ARG}" && -z "${NPS_ANGLE_DEG_ARG}" ) ]]; then
  if [[ ! -f "${KINEMATICS_CONFIG}" ]]; then
    echo "[ERROR] Simulation kinematics config not found: ${KINEMATICS_CONFIG}" >&2
    exit 1
  fi

  KIN_CONFIG_KEY="${KIN}"
  if [[ "${KIN_CONFIG_KEY}" != KinC_* ]]; then
    KIN_CONFIG_KEY="KinC_${KIN_CONFIG_KEY}"
  fi
  KINEMATICS_ROW="$(awk -F',' -v wanted="${KIN_CONFIG_KEY}" '
    function trim(s) { gsub(/^[[:space:]]+|[[:space:]]+$/, "", s); return s }
    NR == 1 {
      for (i = 1; i <= NF; ++i) {
        name = trim($i)
        if (name == "kin_old") kin_col = i
        else if (name == "ebeam_gev") beam_col = i
        else if (name == "nps_theta_deg") theta_col = i
        else if (name == "nps_target_distance_cm") z_col = i
      }
      next
    }
    kin_col > 0 && beam_col > 0 && theta_col > 0 && z_col > 0 && trim($kin_col) == wanted {
      print trim($beam_col) "|" trim($z_col) "|" trim($theta_col)
      exit
    }
  ' "${KINEMATICS_CONFIG}")"
  if [[ -z "${KINEMATICS_ROW}" ]]; then
    echo "[ERROR] No simulation kinematics for ${KIN_CONFIG_KEY} in ${KINEMATICS_CONFIG}" >&2
    exit 1
  fi
  IFS='|' read -r CONFIG_BEAM_ENERGY_GEV CONFIG_Z_NPS_CM CONFIG_NPS_THETA_DEG <<< "${KINEMATICS_ROW}"
  if [[ -z "${BEAM_ENERGY_GEV_ARG}" ]]; then
    BEAM_ENERGY_GEV_ARG="${CONFIG_BEAM_ENERGY_GEV}"
  fi
  if [[ -z "${NPS_Z_NPS_CM_ARG}" ]]; then
    NPS_Z_NPS_CM_ARG="${CONFIG_Z_NPS_CM}"
  fi
  if [[ -z "${NPS_THETA_NPS_DEG_ARG}" && -z "${NPS_ANGLE_DEG_ARG}" ]]; then
    NPS_ANGLE_DEG_ARG="${CONFIG_NPS_THETA_DEG}"
  fi
  echo "[config] ${KIN_CONFIG_KEY}: ebeam=${BEAM_ENERGY_GEV_ARG} GeV, z_nps=${NPS_Z_NPS_CM_ARG} cm, physical theta_nps=${CONFIG_NPS_THETA_DEG} deg"
fi

is_float() {
  [[ "$1" =~ ^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$ ]]
}

if [[ -z "${BEAM_ENERGY_GEV_ARG}" ]]; then
  echo "[ERROR] Beam energy was not resolved from the kinematics config or an override." >&2
  exit 1
fi
if ! is_float "${BEAM_ENERGY_GEV_ARG}"; then
  echo "[ERROR] --beam-energy-gev must be numeric (got: ${BEAM_ENERGY_GEV_ARG})" >&2
  exit 1
fi
if [[ -n "${DEAD_BLOCK_RUN}" && ! "${DEAD_BLOCK_RUN}" =~ ^[1-9][0-9]*$ ]]; then
  echo "[ERROR] --dead-block-run must be a positive integer (got: ${DEAD_BLOCK_RUN})" >&2
  exit 1
fi
if [[ -z "${NPS_Z_NPS_CM_ARG}" ]]; then
  echo "[ERROR] NPS z distance was not resolved from the kinematics config or an override." >&2
  exit 1
fi
if ! is_float "${NPS_Z_NPS_CM_ARG}"; then
  echo "[ERROR] --z-nps-cm must be numeric (got: ${NPS_Z_NPS_CM_ARG})" >&2
  exit 1
fi
if [[ -n "${NPS_THETA_NPS_DEG_ARG}" && -n "${NPS_ANGLE_DEG_ARG}" ]]; then
  echo "[ERROR] Use either --theta-nps-deg/NPS_THETA_NPS_DEG or --nps-angle-deg/NPS_NPS_THETA_DEG, not both." >&2
  exit 1
fi
if [[ -z "${NPS_THETA_NPS_DEG_ARG}" && -z "${NPS_ANGLE_DEG_ARG}" ]]; then
  echo "[ERROR] NPS angle was not resolved from the kinematics config or an override." >&2
  exit 1
fi
if [[ -n "${NPS_THETA_NPS_DEG_ARG}" ]]; then
  if ! is_float "${NPS_THETA_NPS_DEG_ARG}"; then
    echo "[ERROR] --theta-nps-deg must be numeric (got: ${NPS_THETA_NPS_DEG_ARG})" >&2
    exit 1
  fi
fi
if [[ -n "${NPS_ANGLE_DEG_ARG}" ]]; then
  if ! is_float "${NPS_ANGLE_DEG_ARG}"; then
    echo "[ERROR] --nps-angle-deg must be numeric (got: ${NPS_ANGLE_DEG_ARG})" >&2
    exit 1
  fi
fi

SMEAR_GEOMETRY_ARGS=(--z-nps-cm "${NPS_Z_NPS_CM_ARG}")
if [[ -n "${NPS_THETA_NPS_DEG_ARG}" ]]; then
  SMEAR_GEOMETRY_ARGS+=(--theta-nps-deg "${NPS_THETA_NPS_DEG_ARG}")
else
  SMEAR_GEOMETRY_ARGS+=(--nps-angle-deg "${NPS_ANGLE_DEG_ARG}")
fi

if ! command -v "${ROOT_CMD}" >/dev/null 2>&1; then
  echo "[ERROR] ROOT executable not found in PATH (${ROOT_CMD})" >&2
  exit 1
fi
if ! command -v root-config >/dev/null 2>&1; then
  echo "[ERROR] root-config not found in PATH" >&2
  exit 1
fi
if ! command -v "${CXX_CMD}" >/dev/null 2>&1; then
  echo "[ERROR] C++ compiler not found in PATH (${CXX_CMD})" >&2
  exit 1
fi
if ! command -v python3 >/dev/null 2>&1 || [[ ! -f "${PUBLISH_HELPER}" ]]; then
  echo "[ERROR] python3 and ${PUBLISH_HELPER} are required for provenance/publication." >&2
  exit 1
fi

if [[ ! -f "${SMEAR_SRC}" ]]; then
  echo "[ERROR] Missing source file: ${SMEAR_SRC}" >&2
  exit 1
fi
if [[ ! -f "${SIMC_MACRO}" ]]; then
  echo "[ERROR] Missing ROOT macro: ${SIMC_MACRO}" >&2
  exit 1
fi
if [[ ! -f "${ACCEPTANCE_CONFIG}" ]]; then
  echo "[ERROR] Acceptance-cuts config not found: ${ACCEPTANCE_CONFIG}" >&2
  exit 1
fi
if [[ -n "${DEAD_BLOCK_RUN}" && ! -f "${DEAD_BLOCK_CONFIG}" ]]; then
  echo "[ERROR] Dead-block config not found: ${DEAD_BLOCK_CONFIG}" >&2
  exit 1
fi
if [[ ! -f "${COMBINED_FILE}" ]]; then
  echo "[ERROR] Missing combined ROOT file: ${COMBINED_FILE}" >&2
  echo "        Run the combine stage first, or pass --combined-file." >&2
  exit 1
fi

if [[ -n "${SIM_INPUT_EXCLUSIVE}" && ! -f "${SIM_INPUT_EXCLUSIVE}" ]]; then
  echo "[ERROR] Exclusive simulation input not found: ${SIM_INPUT_EXCLUSIVE}" >&2
  exit 1
fi
if [[ -n "${SIM_INPUT_SIDIS}" && ! -f "${SIM_INPUT_SIDIS}" ]]; then
  echo "[ERROR] SIDIS simulation input not found: ${SIM_INPUT_SIDIS}" >&2
  exit 1
fi
if [[ -n "${SIM_INPUT_DELTA}" && ! -f "${SIM_INPUT_DELTA}" ]]; then
  echo "[ERROR] Delta simulation input not found: ${SIM_INPUT_DELTA}" >&2
  exit 1
fi
if [[ ! -f "${SIM_HIST_EXCLUSIVE}" &&
      ( -z "${NPS_SIMC_EXCLUSIVE_NORMFAC:-}" || -z "${NPS_SIMC_EXCLUSIVE_NGEN:-}" ) ]]; then
  echo "[ERROR] Exclusive SIMC hist not found: ${SIM_HIST_EXCLUSIVE}" >&2
  exit 1
fi
if [[ ! -f "${SIM_HIST_SIDIS}" &&
      ( -z "${NPS_SIMC_SIDIS_NORMFAC:-}" || -z "${NPS_SIMC_SIDIS_NGEN:-}" ) ]]; then
  echo "[ERROR] SIDIS SIMC hist not found: ${SIM_HIST_SIDIS}" >&2
  exit 1
fi
if [[ ! -f "${SIM_HIST_DELTA}" &&
      ( -z "${NPS_SIMC_DELTA_NORMFAC:-}" || -z "${NPS_SIMC_DELTA_NGEN:-}" ) ]]; then
  echo "[ERROR] Delta SIMC hist not found: ${SIM_HIST_DELTA}" >&2
  exit 1
fi

for intval in NX NY NSMEAR SMEAR_RANDOM_SEED; do
  if ! [[ "${!intval}" =~ ^[0-9]+$ ]]; then
    echo "[ERROR] ${intval} must be a non-negative integer (got: ${!intval})" >&2
    exit 1
  fi
done

mkdir -p "${ROOT_DIR}" "$(dirname "${OUT_FILE}")" "$(dirname "${SECTION_MAP_FILE}")" "$(dirname "${INTERP_FILE}")" "$(dirname "${INPUT_PREVIEW_PDF}")"
if [[ ! -w "${ROOT_DIR}" ]]; then
  echo "[ERROR] Output directory is not writable: ${ROOT_DIR}" >&2
  exit 1
fi

file_identity() {
  local path="$1"
  local resolved content_hash
  resolved="$(realpath -- "${path}")" || return 1
  content_hash="$(sha256sum -- "${path}" | awk '{print $1}')" || return 1
  [[ "${content_hash}" =~ ^[[:xdigit:]]{64}$ ]] || return 1
  printf '%s|sha256:%s' "${resolved}" "${content_hash}"
}

optional_file_identity() {
  local path="$1"
  if [[ -f "${path}" ]]; then
    file_identity "${path}"
  else
    printf '%s' "not-read:${path}"
  fi
}

if ! command -v sha256sum >/dev/null 2>&1; then
  echo "[ERROR] sha256sum is required for current-run provenance." >&2
  exit 1
fi

INPUT_ID_EXCLUSIVE="$(file_identity "${SIM_INPUT_EXCLUSIVE}")"
INPUT_ID_SIDIS="$(file_identity "${SIM_INPUT_SIDIS}")"
INPUT_ID_DELTA="$(file_identity "${SIM_INPUT_DELTA}")"
HIST_ID_EXCLUSIVE="$(optional_file_identity "${SIM_HIST_EXCLUSIVE}")"
HIST_ID_SIDIS="$(optional_file_identity "${SIM_HIST_SIDIS}")"
HIST_ID_DELTA="$(optional_file_identity "${SIM_HIST_DELTA}")"
ACCEPTANCE_CONFIG_IDENTITY="$(file_identity "${ACCEPTANCE_CONFIG}")"
KINEMATICS_CONFIG_IDENTITY="$(optional_file_identity "${KINEMATICS_CONFIG}")"
DEAD_BLOCK_CONFIG_IDENTITY="$(optional_file_identity "${DEAD_BLOCK_CONFIG}")"
PRODUCER_SOURCE_IDENTITY="pending-compile-dependencies"

input_snapshot() {
  printf '%s\n' \
    "kin=${KIN}" \
    "mode_hint=${MODE_HINT}" \
    "beam_energy_gev=${BEAM_ENERGY_GEV_ARG}" \
    "z_nps_cm=${NPS_Z_NPS_CM_ARG}" \
    "theta_nps_signed_deg=${NPS_THETA_NPS_DEG_ARG:-physical:${NPS_ANGLE_DEG_ARG}}" \
    "dead_block_run=${DEAD_BLOCK_RUN:-inactive}" \
    "requested_smear_seed=${SMEAR_RANDOM_SEED}" \
    "exclusive=${INPUT_ID_EXCLUSIVE}" \
    "sidis=${INPUT_ID_SIDIS}" \
    "delta_pi0=${INPUT_ID_DELTA}" \
    "exclusive_hist=${HIST_ID_EXCLUSIVE}" \
    "sidis_hist=${HIST_ID_SIDIS}" \
    "delta_hist=${HIST_ID_DELTA}" \
    "exclusive_normfac=${NPS_SIMC_EXCLUSIVE_NORMFAC:-from_hist}" \
    "exclusive_ngen=${NPS_SIMC_EXCLUSIVE_NGEN:-from_hist}" \
    "sidis_normfac=${NPS_SIMC_SIDIS_NORMFAC:-from_hist}" \
    "sidis_ngen=${NPS_SIMC_SIDIS_NGEN:-from_hist}" \
    "delta_normfac=${NPS_SIMC_DELTA_NORMFAC:-from_hist}" \
    "delta_ngen=${NPS_SIMC_DELTA_NGEN:-from_hist}" \
    "acceptance=${ACCEPTANCE_CONFIG_IDENTITY}" \
    "kinematics=${KINEMATICS_CONFIG_IDENTITY}" \
    "dead_block=${DEAD_BLOCK_CONFIG_IDENTITY}" \
    "producer_source=${PRODUCER_SOURCE_IDENTITY}"
}

current_input_snapshot() {
  local INPUT_ID_EXCLUSIVE INPUT_ID_SIDIS INPUT_ID_DELTA
  local HIST_ID_EXCLUSIVE HIST_ID_SIDIS HIST_ID_DELTA
  local ACCEPTANCE_CONFIG_IDENTITY KINEMATICS_CONFIG_IDENTITY DEAD_BLOCK_CONFIG_IDENTITY
  local PRODUCER_SOURCE_IDENTITY
  INPUT_ID_EXCLUSIVE="$(file_identity "${SIM_INPUT_EXCLUSIVE}")" || return 1
  INPUT_ID_SIDIS="$(file_identity "${SIM_INPUT_SIDIS}")" || return 1
  INPUT_ID_DELTA="$(file_identity "${SIM_INPUT_DELTA}")" || return 1
  HIST_ID_EXCLUSIVE="$(optional_file_identity "${SIM_HIST_EXCLUSIVE}")" || return 1
  HIST_ID_SIDIS="$(optional_file_identity "${SIM_HIST_SIDIS}")" || return 1
  HIST_ID_DELTA="$(optional_file_identity "${SIM_HIST_DELTA}")" || return 1
  ACCEPTANCE_CONFIG_IDENTITY="$(file_identity "${ACCEPTANCE_CONFIG}")" || return 1
  KINEMATICS_CONFIG_IDENTITY="$(optional_file_identity "${KINEMATICS_CONFIG}")" || return 1
  DEAD_BLOCK_CONFIG_IDENTITY="$(optional_file_identity "${DEAD_BLOCK_CONFIG}")" || return 1
  PRODUCER_SOURCE_IDENTITY="$(python3 "${PUBLISH_HELPER}" --dependency-identity \
    --depfile "${BUILD_DIR}/producer.d" --extra-source "${BASH_SOURCE[0]}" \
    --extra-source "${PUBLISH_HELPER}")" || return 1
  input_snapshot
}

RUN_WORK_DIR="$(mktemp -d "${ROOT_DIR}/.smearing_pipeline.${RUN_TAG}.XXXXXX")"
BUILD_DIR="${RUN_WORK_DIR}/build"
mkdir -p "${BUILD_DIR}"
SMEAR_BIN="${BUILD_DIR}/nps_sim_smearing_new"
SIMC_BIN="${BUILD_DIR}/simc_pi0_analysis"
RUN_SIM_FILE="${RUN_WORK_DIR}/simc_pi0_analysis_output.root"
RUN_SIM_SMEARED_FILE="${RUN_WORK_DIR}/simc_pi0_analysis_output_smeared.root"
RUN_OUT_FILE="${RUN_WORK_DIR}/$(basename "${OUT_FILE}")"
RUN_SECTION_MAP_FILE="${RUN_WORK_DIR}/section_map.csv"
RUN_INTERP_FILE="$(derive_interp_path "${RUN_OUT_FILE}")"
RUN_INPUT_PREVIEW_PDF="${RUN_WORK_DIR}/$(basename "${INPUT_PREVIEW_PDF}")"
ARCHIVE_TMP=""

cleanup_pipeline_staging() {
  if [[ -n "${ARCHIVE_TMP}" && -d "${ARCHIVE_TMP}" ]]; then
    rm -rf -- "${ARCHIVE_TMP}"
  fi
  if [[ -d "${RUN_WORK_DIR}" ]]; then
    rm -rf -- "${RUN_WORK_DIR}"
  fi
}
trap cleanup_pipeline_staging EXIT

# Freeze small inputs for both producer passes. Original paths/content hashes
# remain in the provenance snapshot and are checked again before publication.
mkdir -p "${RUN_WORK_DIR}/inputs"
freeze_optional_input() {
  local original="$1" name="$2" expected="$3"
  if [[ ! -f "${original}" ]]; then
    printf '%s' "${original}"
    return
  fi
  local frozen="${RUN_WORK_DIR}/inputs/${name}"
  cp -p -- "${original}" "${frozen}"
  if [[ "$(sha256sum "${frozen}" | awk '{print $1}')" != "${expected##*|sha256:}" ]]; then
    echo "[ERROR] Input changed while snapshotting: ${original}" >&2
    return 1
  fi
  printf '%s' "${frozen}"
}
FROZEN_ACCEPTANCE_CONFIG="$(freeze_optional_input "${ACCEPTANCE_CONFIG}" acceptance_cuts.conf "${ACCEPTANCE_CONFIG_IDENTITY}")"
FROZEN_KINEMATICS_CONFIG="$(freeze_optional_input "${KINEMATICS_CONFIG}" kinematics.csv "${KINEMATICS_CONFIG_IDENTITY}")"
FROZEN_DEAD_BLOCK_CONFIG="$(freeze_optional_input "${DEAD_BLOCK_CONFIG}" dead_blocks.csv "${DEAD_BLOCK_CONFIG_IDENTITY}")"
FROZEN_HIST_EXCLUSIVE="$(freeze_optional_input "${SIM_HIST_EXCLUSIVE}" exclusive.hist "${HIST_ID_EXCLUSIVE}")"
FROZEN_HIST_SIDIS="$(freeze_optional_input "${SIM_HIST_SIDIS}" sidis.hist "${HIST_ID_SIDIS}")"
FROZEN_HIST_DELTA="$(freeze_optional_input "${SIM_HIST_DELTA}" delta.hist "${HIST_ID_DELTA}")"

run_simc_analysis() {
  local stage_label="$1"
  local smearing_mode="$2"
  local producer_mode="$3"
  local nominal_identity="${4:-}"
  local -a env_args=(
    "NPS_KIN=${KIN}"
    "NPS_SMEAR_KIN=${KIN}"
    "NPS_MODE=${MODE_HINT}"
    "NPS_SMEAR_MODE_HINT=${MODE_HINT}"
    "NPS_ACCEPTANCE_CUTS_CONFIG=${FROZEN_ACCEPTANCE_CONFIG}"
    "NPS_SIMULATION_KINEMATICS_CONFIG=${FROZEN_KINEMATICS_CONFIG}"
    "NPS_SIMC_PRODUCTION_DIR=${SIMC_PRODUCTION_DIR}"
    "NPS_GEANT_PRODUCTION_DIR=${GEANT_PRODUCTION_DIR}"
    "NPS_DEAD_BLOCK_CONFIG=${FROZEN_DEAD_BLOCK_CONFIG}"
    "NPS_SMEAR_FILE=${RUN_OUT_FILE}"
    "NPS_SMEAR_INTERP_FILE=${RUN_INTERP_FILE}"
    "NPS_SECTION_MAP_FILE=${RUN_SECTION_MAP_FILE}"
    "NPS_Z_NPS_CM=${NPS_Z_NPS_CM_ARG}"
    "NPS_SIMC_OUTPUT_FILE=${RUN_SIM_FILE}"
    "NPS_SIMC_OUTPUT_FILE_SMEARED=${RUN_SIM_SMEARED_FILE}"
    "NPS_SIMC_SMEARED_OUTPUT_FILE=${RUN_SIM_SMEARED_FILE}"
    "NPS_SIMC_PRODUCTION_MODE=${producer_mode}"
    "NPS_SMEAR_RANDOM_SEED=${SMEAR_RANDOM_SEED}"
    "NPS_EBEAM=${BEAM_ENERGY_GEV_ARG}"
    "NPS_PIPELINE_RUN_ID=${RUN_TAG}"
    "NPS_PRODUCER_INPUT_SET_IDENTITY=${PRODUCER_INPUT_SET_IDENTITY}"
    "NPS_PRODUCER_SOURCE_IDENTITY=${PRODUCER_SOURCE_IDENTITY}"
    "NPS_ACCEPTANCE_CONFIG_IDENTITY=${ACCEPTANCE_CONFIG_IDENTITY}"
    "NPS_KINEMATICS_CONFIG_IDENTITY=${KINEMATICS_CONFIG_IDENTITY}"
    "NPS_DEAD_BLOCK_CONFIG_IDENTITY=${DEAD_BLOCK_CONFIG_IDENTITY}"
    "NPS_INPUT_IDENTITY_exclusive=${INPUT_ID_EXCLUSIVE}"
    "NPS_INPUT_IDENTITY_sidis=${INPUT_ID_SIDIS}"
    "NPS_INPUT_IDENTITY_delta_pi0=${INPUT_ID_DELTA}"
    "NPS_NOMINAL_INPUT_IDENTITY=${nominal_identity}"
  )

  if [[ -n "${NPS_THETA_NPS_DEG_ARG}" ]]; then
    env_args+=("NPS_THETA_NPS_DEG=${NPS_THETA_NPS_DEG_ARG}")
  else
    env_args+=("NPS_NPS_THETA_DEG=${NPS_ANGLE_DEG_ARG}")
  fi

  if [[ -n "${smearing_mode}" ]]; then
    env_args+=("NPS_SMEARING_MODE=${smearing_mode}")
  fi

  env_args+=(
    "NPS_SIMC_EXCLUSIVE_INPUT=${SIM_INPUT_EXCLUSIVE}"
    "NPS_SIMC_INPUT_FILE_EXCLUSIVE=${SIM_INPUT_EXCLUSIVE}"
    "NPS_SIMC_SIDIS_INPUT=${SIM_INPUT_SIDIS}"
    "NPS_SIMC_INPUT_FILE_SIDIS=${SIM_INPUT_SIDIS}"
    "NPS_SIMC_DELTA_INPUT=${SIM_INPUT_DELTA}"
    "NPS_SIMC_INPUT_FILE_DELTA=${SIM_INPUT_DELTA}"
    "NPS_SIMC_EXCLUSIVE_HIST=${FROZEN_HIST_EXCLUSIVE}"
    "NPS_SIMC_SIDIS_HIST=${FROZEN_HIST_SIDIS}"
    "NPS_SIMC_DELTA_HIST=${FROZEN_HIST_DELTA}"
  )
  if [[ -n "${DEAD_BLOCK_RUN}" ]]; then
    env_args+=("NPS_SIMC_DEAD_BLOCK_RUN=${DEAD_BLOCK_RUN}")
  fi

  echo ""
  echo "============================================================================"
  echo "${stage_label}"
  echo "============================================================================"
  set +e
  env "${env_args[@]}" "${SIMC_BIN}"
  local status=$?
  set -e
  if [[ ${status} -ne 0 ]]; then
    echo "[ERROR] SIMC producer failed in ${producer_mode} mode with status ${status}." >&2
  fi
  return "${status}"
}

echo ""
echo "============================================================================"
echo "Step 1/5: Compiling producer and fitter"
echo "============================================================================"
"${CXX_CMD}" "${SIMC_MACRO}" $(root-config --cflags) -O2 -std=c++17 -I"${REPO_ROOT}/src" -MM -MF "${BUILD_DIR}/producer.d"
"${CXX_CMD}" "${SMEAR_SRC}" $(root-config --cflags) -O3 -march=native -std=c++17 -fopenmp -I"${REPO_ROOT}/src" -MM -MF "${BUILD_DIR}/fitter.d"
PRODUCER_SOURCE_BEFORE="$(python3 "${PUBLISH_HELPER}" --dependency-identity \
  --depfile "${BUILD_DIR}/producer.d" --extra-source "${BASH_SOURCE[0]}" \
  --extra-source "${PUBLISH_HELPER}")"
FITTER_SOURCE_BEFORE="$(python3 "${PUBLISH_HELPER}" --dependency-identity --depfile "${BUILD_DIR}/fitter.d")"
PROTECTED_INPUT_ARGS=()
for protected in "${SIM_INPUT_EXCLUSIVE}" "${SIM_INPUT_SIDIS}" "${SIM_INPUT_DELTA}" \
                 "${SIM_HIST_EXCLUSIVE}" "${SIM_HIST_SIDIS}" "${SIM_HIST_DELTA}" \
                 "${COMBINED_FILE}" "${ACCEPTANCE_CONFIG}" "${KINEMATICS_CONFIG}" \
                 "${DEAD_BLOCK_CONFIG}"; do
  PROTECTED_INPUT_ARGS+=(--protected-path "${protected}")
done
LATEST_MANIFEST="$(dirname "${OUT_FILE}")/smearing_latest.json"
python3 "${PUBLISH_HELPER}" --check-destinations \
  --destination "${SIM_FILE}" --destination "${SIM_SMEARED_FILE}" \
  --destination "${OUT_FILE}" --destination "${SECTION_MAP_FILE}" \
  --destination "${INTERP_FILE}" --destination "${INPUT_PREVIEW_PDF}" \
  --destination "${LATEST_MANIFEST}" "${PROTECTED_INPUT_ARGS[@]}" \
  --depfile "${BUILD_DIR}/producer.d" --depfile "${BUILD_DIR}/fitter.d" \
  --extra-source "${BASH_SOURCE[0]}" --extra-source "${PUBLISH_HELPER}"
"${CXX_CMD}" "${SIMC_MACRO}" $(root-config --cflags --libs) -O2 -std=c++17 -I"${REPO_ROOT}/src" -MMD -MF "${BUILD_DIR}/producer.d" -o "${SIMC_BIN}"
"${CXX_CMD}" "${SMEAR_SRC}" $(root-config --cflags --libs) -lMathMore -O3 -march=native -std=c++17 -fopenmp -I"${REPO_ROOT}/src" -MMD -MF "${BUILD_DIR}/fitter.d" -o "${SMEAR_BIN}"
PRODUCER_SOURCE_IDENTITY="$(python3 "${PUBLISH_HELPER}" --dependency-identity \
  --depfile "${BUILD_DIR}/producer.d" --extra-source "${BASH_SOURCE[0]}" \
  --extra-source "${PUBLISH_HELPER}")"
if [[ "${PRODUCER_SOURCE_IDENTITY}" != "${PRODUCER_SOURCE_BEFORE}" || \
      "$(python3 "${PUBLISH_HELPER}" --dependency-identity --depfile "${BUILD_DIR}/fitter.d")" != "${FITTER_SOURCE_BEFORE}" ]]; then
  echo "[ERROR] Source/dependency content changed during compilation; refusing mismatched build provenance." >&2
  exit 1
fi
{
  "${CXX_CMD}" --version
  root-config --version --cflags --libs
  printf 'producer_flags=-O2 -std=c++17 -I%s/src\n' "${REPO_ROOT}"
  printf 'fitter_flags=-lMathMore -O3 -march=native -std=c++17 -fopenmp -I%s/src\n' "${REPO_ROOT}"
} > "${BUILD_DIR}/toolchain.txt"
python3 "${PUBLISH_HELPER}" --build-provenance "${RUN_WORK_DIR}/build_provenance.json" \
  --depfile "${BUILD_DIR}/producer.d" --depfile "${BUILD_DIR}/fitter.d" \
  --extra-source "${BASH_SOURCE[0]}" --extra-source "${PUBLISH_HELPER}" \
  --executable "${SIMC_BIN}" --executable "${SMEAR_BIN}" \
  --toolchain-info "${BUILD_DIR}/toolchain.txt"
INPUT_SNAPSHOT_BEFORE="$(input_snapshot)"
PRODUCER_INPUT_SET_IDENTITY="$(printf '%s' "${INPUT_SNAPSHOT_BEFORE}" | sha256sum | awk '{print $1}')"
if [[ "$(current_input_snapshot)" != "${INPUT_SNAPSHOT_BEFORE}" ]]; then
  echo "[ERROR] Producer inputs/configuration changed during preparation." >&2
  exit 1
fi

run_simc_analysis "Step 2/5: Generating fresh nominal simulation" "off" "nominal-only"
"${SIMC_BIN}" --validate-output "${RUN_SIM_FILE}" nominal "${RUN_TAG}" \
  "${PRODUCER_INPUT_SET_IDENTITY}" "" -1
NOMINAL_ENTRY_COUNT="$("${SIMC_BIN}" --entry-count "${RUN_SIM_FILE}")"
NOMINAL_INPUT_IDENTITY="$(sha256sum "${RUN_SIM_FILE}" | awk '{print $1}')"

echo ""
echo "============================================================================"
echo "Step 3/5: Running section smearing fit"
echo "============================================================================"
echo "The fitter will first write all-observable input overlays to:"
echo "  ${RUN_INPUT_PREVIEW_PDF}"
echo "Optimization will remain paused until those plots are explicitly approved."
set +e
NPS_PIPELINE_RUN_ID="${RUN_TAG}" NPS_NOMINAL_INPUT_IDENTITY="${NOMINAL_INPUT_IDENTITY}" "${SMEAR_BIN}" \
  "${COMBINED_FILE}" "${DATA_TREE}" \
  "${RUN_SIM_FILE}" "${SIM_TREE}" \
  "${RUN_OUT_FILE}" "${NX}" "${NY}" "${X_MIN}" "${X_MAX}" "${Y_MIN}" "${Y_MAX}" "${OVERLAP}" "${NSMEAR}" \
  --beam-energy-gev "${BEAM_ENERGY_GEV_ARG}" \
  "${SMEAR_GEOMETRY_ARGS[@]}" \
  --input-preview-pdf "${RUN_INPUT_PREVIEW_PDF}" \
  --confirm-inputs
SMEAR_STATUS=$?
set -e
if [[ ${SMEAR_STATUS} -eq 20 ]]; then
  if [[ -s "${RUN_INPUT_PREVIEW_PDF}" ]]; then
    CANCEL_PREVIEW_TMP="$(mktemp "${INPUT_PREVIEW_PDF}.tmp.XXXXXX")"
    cp -p -- "${RUN_INPUT_PREVIEW_PDF}" "${CANCEL_PREVIEW_TMP}"
    mv -f -- "${CANCEL_PREVIEW_TMP}" "${INPUT_PREVIEW_PDF}"
    echo "Input preview retained at: ${INPUT_PREVIEW_PDF}"
  fi
  echo "Pipeline stopped before optimization because the input histograms were not approved."
  exit 0
fi
if [[ ${SMEAR_STATUS} -ne 0 ]]; then
  echo "[ERROR] Smearing fitter exited with status ${SMEAR_STATUS}." >&2
  exit "${SMEAR_STATUS}"
fi

if [[ ! -s "${RUN_OUT_FILE}" ]]; then
  echo "[ERROR] Fitter ROOT output was not generated for this run: ${RUN_OUT_FILE}" >&2
  exit 1
fi
if [[ ! -s "${RUN_SECTION_MAP_FILE}" ]]; then
  echo "[ERROR] Section map was not generated for this run: ${RUN_SECTION_MAP_FILE}" >&2
  exit 1
fi
if [[ ! -s "${RUN_INTERP_FILE}" ]]; then
  echo "[ERROR] Interpolated map ROOT was not generated for this run: ${RUN_INTERP_FILE}" >&2
  exit 1
fi
RUN_LOOKUP_COMPARISON_CSV="${RUN_WORK_DIR}/smearing_response_lookup_comparison.csv"
if ! awk -F',' '
  NR == 1 {ok = ($1 == "observable" && $2 == "n_events" && $3 == "nsmear")}
  NR > 1 && ($1 == "mgg" || $1 == "mmiss" || $1 == "mpgg2") {seen[$1] = 1}
  END {exit !(ok && seen["mgg"] && seen["mmiss"] && seen["mpgg2"])}
' "${RUN_LOOKUP_COMPARISON_CSV}"; then
  echo "[ERROR] Current-run response-lookup comparison is missing or invalid: ${RUN_LOOKUP_COMPARISON_CSV}" >&2
  exit 1
fi
if ! awk -F',' 'NR == 1 {ok = ($1 == "ix" && $2 == "iy")} NR > 1 && ($NF == "fit_ok" || $NF == "poor_fit") {rows++} END {exit !(ok && rows > 0)}' "${RUN_SECTION_MAP_FILE}"; then
  echo "[ERROR] Current-run section map has no usable fitted rows: ${RUN_SECTION_MAP_FILE}" >&2
  exit 1
fi

env NPS_VALIDATE_FIT_FILE="${RUN_OUT_FILE}" \
    NPS_VALIDATE_INTERP_FILE="${RUN_INTERP_FILE}" \
    NPS_VALIDATE_RUN_ID="${RUN_TAG}" \
    NPS_VALIDATE_NOMINAL_ID="${NOMINAL_INPUT_IDENTITY}" \
    "${ROOT_CMD}" -l -b -q -e '
      int rc = 0;
      TFile fit(gSystem->Getenv("NPS_VALIDATE_FIT_FILE"), "READ");
      TFile interp(gSystem->Getenv("NPS_VALIDATE_INTERP_FILE"), "READ");
      if (fit.IsZombie() || fit.TestBit(TFile::kRecovered)) rc = 21;
      TNamed *run = rc ? nullptr : dynamic_cast<TNamed*>(fit.Get("run_tag"));
      if (!rc && (!run || TString(run->GetTitle()) != gSystem->Getenv("NPS_VALIDATE_RUN_ID"))) rc = 22;
      TNamed *nominal = rc ? nullptr : dynamic_cast<TNamed*>(fit.Get("nominal_input_identity"));
      if (!rc && (!nominal || TString(nominal->GetTitle()) != gSystem->Getenv("NPS_VALIDATE_NOMINAL_ID"))) rc = 22;
      const char *maps[] = {"h_mu_a_interp", "h_mu_b_interp", "h_mu_c_interp", "h_sigma_interp", "h_sigma_pos_interp"};
      if (!rc && (interp.IsZombie() || interp.TestBit(TFile::kRecovered))) rc = 23;
      for (const char *name : maps) if (!rc && !interp.Get(name)) rc = 24;
      TNamed *interp_run = rc ? nullptr : dynamic_cast<TNamed*>(interp.Get("pipeline_run_id"));
      TNamed *interp_nominal = rc ? nullptr : dynamic_cast<TNamed*>(interp.Get("nominal_input_identity"));
      if (!rc && (!interp_run || TString(interp_run->GetTitle()) != gSystem->Getenv("NPS_VALIDATE_RUN_ID") ||
                  !interp_nominal || TString(interp_nominal->GetTitle()) != gSystem->Getenv("NPS_VALIDATE_NOMINAL_ID"))) rc = 25;
      gSystem->Exit(rc);'

if [[ "$(current_input_snapshot)" != "${INPUT_SNAPSHOT_BEFORE}" ]]; then
  echo "[ERROR] Producer input/configuration identity changed after nominal generation; refusing final production." >&2
  exit 1
fi

SECTION_MAP_IDENTITY="$(sha256sum "${RUN_SECTION_MAP_FILE}" | awk '{print $1}')"
INTERP_MAP_IDENTITY="$(sha256sum "${RUN_INTERP_FILE}" | awk '{print $1}')"
export NPS_SECTION_MAP_IDENTITY="${SECTION_MAP_IDENTITY}"
export NPS_INTERP_MAP_IDENTITY="${INTERP_MAP_IDENTITY}"

run_simc_analysis "Step 4/5: Producing smeared-only simulation from the same raw inputs" \
  "section" "smeared-only" "${NOMINAL_INPUT_IDENTITY}"
"${SIMC_BIN}" --validate-output "${RUN_SIM_SMEARED_FILE}" smeared "${RUN_TAG}" \
  "${PRODUCER_INPUT_SET_IDENTITY}" "${NOMINAL_INPUT_IDENTITY}" "${NOMINAL_ENTRY_COUNT}"

if [[ "$(sha256sum "${RUN_SIM_FILE}" | awk '{print $1}')" != "${NOMINAL_INPUT_IDENTITY}" ]]; then
  echo "[ERROR] Fresh nominal fit input changed during fitting/final production." >&2
  exit 1
fi
if [[ "$(current_input_snapshot)" != "${INPUT_SNAPSHOT_BEFORE}" ]]; then
  echo "[ERROR] Producer input/configuration identity changed during final production." >&2
  exit 1
fi

RUN_ARCHIVE_DIR="$(dirname "${OUT_FILE}")/smearing_runs/${RUN_TAG}"
if [[ -e "${RUN_ARCHIVE_DIR}" ]]; then
  echo "[ERROR] Refusing to replace existing run archive: ${RUN_ARCHIVE_DIR}" >&2
  exit 1
fi
mkdir -p "$(dirname "${RUN_ARCHIVE_DIR}")"
ARCHIVE_TMP="$(mktemp -d "$(dirname "${RUN_ARCHIVE_DIR}")/.${RUN_TAG}.publish.XXXXXX")"

archive_current_file() {
  local src="$1" name="$2"
  [[ -f "${src}" ]] || return 0
  if [[ -e "${ARCHIVE_TMP}/${name}" ]]; then
    echo "[ERROR] Archive artifact name collision: ${name}" >&2
    return 1
  fi
  cp -p -- "${src}" "${ARCHIVE_TMP}/${name}"
}

archive_current_file "${RUN_SIM_FILE}" "simc_pi0_analysis_output.root"
archive_current_file "${RUN_SIM_SMEARED_FILE}" "simc_pi0_analysis_output_smeared.root"
archive_current_file "${RUN_OUT_FILE}" "$(basename "${OUT_FILE}")"
archive_current_file "${RUN_SECTION_MAP_FILE}" "$(basename "${SECTION_MAP_FILE}")"
archive_current_file "${RUN_INTERP_FILE}" "$(basename "${INTERP_FILE}")"
archive_current_file "${RUN_INPUT_PREVIEW_PDF}" "$(basename "${INPUT_PREVIEW_PDF}")"
archive_current_file "${RUN_WORK_DIR}/build_provenance.json" "build_provenance.json"
cp -a -- "${RUN_WORK_DIR}/inputs" "${ARCHIVE_TMP}/inputs"
for name in smearing_optimizer_summary.csv smearing_optimizer_seeds.csv \
            smearing_optimizer_profiles.csv smearing_closure_summary.csv \
            smearing_sweep_history.csv smearing_objective_breakdown.csv \
            smearing_response_lookup_comparison.csv smearing_config_fingerprint.txt \
            chi2_scans.pdf; do
  archive_current_file "${RUN_WORK_DIR}/${name}" "${name}"
done
if [[ -d "${RUN_WORK_DIR}/smearing_runs/${RUN_TAG}/progress" ]]; then
  cp -a -- "${RUN_WORK_DIR}/smearing_runs/${RUN_TAG}/progress" "${ARCHIVE_TMP}/progress"
fi
for ending in .pdf .png .root _params.csv _peak_scan.csv; do
  archive_current_file "${RUN_WORK_DIR}/simc_unsmeared_2d_mass_cut${ending}" \
                       "simc_unsmeared_2d_mass_cut${ending}"
  archive_current_file "${RUN_WORK_DIR}/simc_smeared_2d_mass_cut${ending}" \
                       "simc_smeared_2d_mass_cut${ending}"
done

PROVENANCE_FILE="${ARCHIVE_TMP}/pipeline_provenance.txt"
{
  printf 'pipeline_run_id=%s\n' "${RUN_TAG}"
  printf 'nominal_input_sha256=%s\n' "${NOMINAL_INPUT_IDENTITY}"
  printf 'producer_input_set_sha256=%s\n' "${PRODUCER_INPUT_SET_IDENTITY}"
  printf 'section_map_sha256=%s\n' "${SECTION_MAP_IDENTITY}"
  printf 'interpolated_map_sha256=%s\n' "${INTERP_MAP_IDENTITY}"
  printf 'combined_input=%s\n' "$(file_identity "${COMBINED_FILE}")"
  printf 'beam_energy_gev=%s\n' "${BEAM_ENERGY_GEV_ARG}"
  printf 'nps_z_cm=%s\n' "${NPS_Z_NPS_CM_ARG}"
  printf 'nps_theta_signed_deg=%s\n' "${NPS_THETA_NPS_DEG_ARG:-physical:${NPS_ANGLE_DEG_ARG}}"
  printf 'mode_hint=%s\n' "${MODE_HINT}"
  printf 'dead_block_run=%s\n' "${DEAD_BLOCK_RUN:-inactive}"
  printf '%s\n' "${INPUT_SNAPSHOT_BEFORE}"
} > "${PROVENANCE_FILE}"

PUBLICATION_ARGS=()
queue_publication() {
  local src="$1" dst="$2" archive_name="${3:-$(basename "$1")}"
  PUBLICATION_ARGS+=(--file "${ARCHIVE_TMP}/${archive_name}" "${dst}")
}

echo ""
echo "============================================================================"
echo "Step 5/5: Publishing validated current-run artifacts"
echo "============================================================================"
queue_publication "${RUN_SIM_FILE}" "${SIM_FILE}"
queue_publication "${RUN_SIM_SMEARED_FILE}" "${SIM_SMEARED_FILE}"
queue_publication "${RUN_OUT_FILE}" "${OUT_FILE}"
queue_publication "${RUN_SECTION_MAP_FILE}" "${SECTION_MAP_FILE}" "$(basename "${SECTION_MAP_FILE}")"
queue_publication "${RUN_INTERP_FILE}" "${INTERP_FILE}" "$(basename "${INTERP_FILE}")"
queue_publication "${RUN_INPUT_PREVIEW_PDF}" "${INPUT_PREVIEW_PDF}"
for name in smearing_optimizer_summary.csv smearing_optimizer_seeds.csv \
            smearing_optimizer_profiles.csv smearing_closure_summary.csv \
            smearing_sweep_history.csv smearing_objective_breakdown.csv \
            smearing_response_lookup_comparison.csv smearing_config_fingerprint.txt \
            chi2_scans.pdf; do
  if [[ -f "${RUN_WORK_DIR}/${name}" ]]; then
    queue_publication "${RUN_WORK_DIR}/${name}" "$(dirname "${OUT_FILE}")/${name}"
  fi
done
for ending in .pdf .png .root _params.csv _peak_scan.csv; do
  if [[ -f "${RUN_WORK_DIR}/simc_unsmeared_2d_mass_cut${ending}" ]]; then
    queue_publication "${RUN_WORK_DIR}/simc_unsmeared_2d_mass_cut${ending}" \
                 "$(dirname "${SIM_FILE}")/simc_unsmeared_2d_mass_cut${ending}"
  fi
  if [[ -f "${RUN_WORK_DIR}/simc_smeared_2d_mass_cut${ending}" ]]; then
    queue_publication "${RUN_WORK_DIR}/simc_smeared_2d_mass_cut${ending}" \
                 "$(dirname "${SIM_SMEARED_FILE}")/simc_smeared_2d_mass_cut${ending}"
  fi
done
if [[ "$(sha256sum "${ARCHIVE_TMP}/simc_pi0_analysis_output.root" | awk '{print $1}')" != "${NOMINAL_INPUT_IDENTITY}" ]]; then
  echo "[ERROR] Archived nominal copy is not byte-identical to the fitter input." >&2
  exit 1
fi
python3 "${PUBLISH_HELPER}" --archive-staging "${ARCHIVE_TMP}" \
  --archive "${RUN_ARCHIVE_DIR}" --manifest "${LATEST_MANIFEST}" \
  --run-id "${RUN_TAG}" "${PUBLICATION_ARGS[@]}" "${PROTECTED_INPUT_ARGS[@]}" \
  --depfile "${BUILD_DIR}/producer.d" --depfile "${BUILD_DIR}/fitter.d" \
  --extra-source "${BASH_SOURCE[0]}" --extra-source "${PUBLISH_HELPER}"
ARCHIVE_TMP=""

echo ""
echo "============================================================================"
echo "Pipeline complete"
echo "  kin:             ${KIN}"
echo "  target:          ${TARGET}"
echo "  mode hint:       ${MODE_HINT}"
echo "  requested seed:  ${SMEAR_RANDOM_SEED} (recorded; producer RNG policy unchanged)"
echo "  beam energy GeV: ${BEAM_ENERGY_GEV_ARG}"
echo "  nps z [cm]:      ${NPS_Z_NPS_CM_ARG}"
if [[ -n "${NPS_THETA_NPS_DEG_ARG}" ]]; then
  echo "  nps theta [deg]: ${NPS_THETA_NPS_DEG_ARG} (NPS->Hall)"
else
  echo "  nps angle [deg]: ${NPS_ANGLE_DEG_ARG} (physical; fitter uses negative)"
fi
echo "  acceptance cfg:  ${ACCEPTANCE_CONFIG}"
echo "  combined input:  ${COMBINED_FILE}"
echo "  fresh sim input: ${SIM_FILE} (sha256=${NOMINAL_INPUT_IDENTITY})"
echo "  smear fit root:  ${OUT_FILE}"
echo "  section map:     ${SECTION_MAP_FILE}"
echo "  interp map:      ${INTERP_FILE}"
echo "  input preview:   ${INPUT_PREVIEW_PDF}"
echo "  sim smeared out: ${SIM_SMEARED_FILE}"
echo "  run archive:     ${RUN_ARCHIVE_DIR}"
echo "  latest manifest: ${LATEST_MANIFEST} (use archived paths for a consistent set)"
echo "============================================================================"
