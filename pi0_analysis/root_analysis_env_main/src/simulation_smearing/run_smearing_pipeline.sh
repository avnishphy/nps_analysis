#!/usr/bin/env bash
# ============================================================================
# End-to-end NPS smearing pipeline
# 1) Ensure unsmeared SIMC analysis output exists
# 2) Compile and run nps_sim_smearing_new.C (section fit + interpolated maps)
# 3) Regenerate smeared SIMC output via simc_pi0_analysis.C
# ============================================================================

# ./src/simulation_smearing/run_smearing_pipeline.sh --kin KinC_x36_5 --target LH2 --combined-file output/KinC_x36_5_yaopeng_08gev/KinC_x36_5/root/combined_branches_LH2_minus_3728_3729.root --root-dir output/KinC_x36_5_yaopeng_08gev/KinC_x36_5/root --nx 6 --ny 6 --x-min -27 --x-max 27 --y-min -33 --y-max 33 --overlap 0.1 --beam-energy-gev 10.538 --z-nps-cm 407 --nps-angle-deg 12.2

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

ROOT_CMD="${ROOT_CMD:-root}"
CXX_CMD="${CXX:-g++}"

SMEAR_SRC="${SCRIPT_DIR}/nps_sim_smearing_new.C"
SIMC_MACRO="${SCRIPT_DIR}/simc_pi0_analysis.C"

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
  --smear-seed <int>          RNG seed for simc smearing producer (default: 42)
  --beam-energy-gev <float>   Required beam energy used by the smearing fitter
  --z-nps-cm <float>          Required NPS target-to-calorimeter distance in cm
  --theta-nps-deg <float>     Required signed NPS->Hall rotation angle in degrees
  --nps-angle-deg <float>     Required physical NPS angle in degrees; fitter uses its negative
  --mode <auto|hcana|waveform>
                               Mode hint for timing defaults in SIM producer (default: auto)
  --acceptance-config <path>  Acceptance-cuts config file (default: repo/config/acceptance_cuts.conf)
  --sim-excl-input <path>     Exclusive SIMC input (forwarded to simc_pi0_analysis.C)
  --sim-sidis-input <path>    SIDIS SIMC input (forwarded to simc_pi0_analysis.C)
  --sim-delta-input <path>    Delta pi0 SIMC input (forwarded to simc_pi0_analysis.C)
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

is_float() {
  [[ "$1" =~ ^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$ ]]
}

if [[ -z "${BEAM_ENERGY_GEV_ARG}" ]]; then
  echo "[ERROR] Beam energy is required. Pass --beam-energy-gev <GeV> or set NPS_EBEAM." >&2
  exit 1
fi
if ! is_float "${BEAM_ENERGY_GEV_ARG}"; then
  echo "[ERROR] --beam-energy-gev must be numeric (got: ${BEAM_ENERGY_GEV_ARG})" >&2
  exit 1
fi
if [[ -z "${NPS_Z_NPS_CM_ARG}" ]]; then
  echo "[ERROR] NPS z distance is required. Pass --z-nps-cm <cm> or set NPS_Z_NPS_CM." >&2
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
  echo "[ERROR] NPS angle is required. Pass --theta-nps-deg <signed_NPS_to_Hall_deg> or --nps-angle-deg <physical_deg>." >&2
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

run_simc_analysis() {
  local stage_label="$1"
  local smearing_mode="${2:-}"
  local -a env_args=(
    "NPS_KIN=${KIN}"
    "NPS_SMEAR_KIN=${KIN}"
    "NPS_MODE=${MODE_HINT}"
    "NPS_SMEAR_MODE_HINT=${MODE_HINT}"
    "NPS_ACCEPTANCE_CUTS_CONFIG=${ACCEPTANCE_CONFIG}"
    "NPS_SMEAR_FILE=${OUT_FILE}"
    "NPS_SMEAR_INTERP_FILE=${INTERP_FILE}"
    "NPS_SECTION_MAP_FILE=${SECTION_MAP_FILE}"
    "NPS_Z_NPS_CM=${NPS_Z_NPS_CM_ARG}"
    "NPS_SIMC_OUTPUT_FILE=${SIM_FILE}"
    "NPS_SIMC_OUTPUT_FILE_SMEARED=${SIM_SMEARED_FILE}"
    "NPS_SIMC_SMEARED_OUTPUT_FILE=${SIM_SMEARED_FILE}"
    "NPS_SMEAR_RANDOM_SEED=${SMEAR_RANDOM_SEED}"
    "NPS_EBEAM=${BEAM_ENERGY_GEV_ARG}"
  )

  if [[ -n "${NPS_THETA_NPS_DEG_ARG}" ]]; then
    env_args+=("NPS_THETA_NPS_DEG=${NPS_THETA_NPS_DEG_ARG}")
  else
    env_args+=("NPS_NPS_THETA_DEG=${NPS_ANGLE_DEG_ARG}")
  fi

  if [[ -n "${smearing_mode}" ]]; then
    env_args+=("NPS_SMEARING_MODE=${smearing_mode}")
  fi

  if [[ -n "${SIM_INPUT_EXCLUSIVE}" ]]; then
    env_args+=("NPS_SIMC_EXCLUSIVE_INPUT=${SIM_INPUT_EXCLUSIVE}" "NPS_SIMC_INPUT_FILE_EXCLUSIVE=${SIM_INPUT_EXCLUSIVE}")
  fi
  if [[ -n "${SIM_INPUT_SIDIS}" ]]; then
    env_args+=("NPS_SIMC_SIDIS_INPUT=${SIM_INPUT_SIDIS}" "NPS_SIMC_INPUT_FILE_SIDIS=${SIM_INPUT_SIDIS}")
  fi
  if [[ -n "${SIM_INPUT_DELTA}" ]]; then
    env_args+=("NPS_SIMC_DELTA_INPUT=${SIM_INPUT_DELTA}" "NPS_SIMC_INPUT_FILE_DELTA=${SIM_INPUT_DELTA}")
  fi

  echo ""
  echo "============================================================================"
  echo "${stage_label}"
  echo "============================================================================"
  env "${env_args[@]}" "${ROOT_CMD}" -l -b -q "${SIMC_MACRO}()"
}

if [[ ! -f "${SIM_FILE}" ]]; then
  run_simc_analysis "Step 1/4: Unsmeared simulation file missing; generating simc output" "off"
  if [[ ! -f "${SIM_FILE}" ]]; then
    echo "[ERROR] Failed to create unsmeared simulation output: ${SIM_FILE}" >&2
    exit 1
  fi
fi

BUILD_DIR="$(mktemp -d "${ROOT_DIR}/.build_smearing.XXXXXX")"
trap 'rm -rf "${BUILD_DIR}"' EXIT
SMEAR_BIN="${BUILD_DIR}/nps_sim_smearing_new"

echo ""
echo "============================================================================"
echo "Step 2/4: Compiling nps_sim_smearing_new.C"
echo "============================================================================"
"${CXX_CMD}" "${SMEAR_SRC}" $(root-config --cflags --libs) -lMathMore -O2 -std=c++17 -fopenmp -I"${REPO_ROOT}/src" -o "${SMEAR_BIN}"

echo ""
echo "============================================================================"
echo "Step 3/4: Running section smearing fit"
echo "============================================================================"
echo "The fitter will first write all-observable input overlays to:"
echo "  ${INPUT_PREVIEW_PDF}"
echo "Optimization will remain paused until those plots are explicitly approved."
set +e
"${SMEAR_BIN}" \
  "${COMBINED_FILE}" "${DATA_TREE}" \
  "${SIM_FILE}" "${SIM_TREE}" \
  "${OUT_FILE}" "${NX}" "${NY}" "${X_MIN}" "${X_MAX}" "${Y_MIN}" "${Y_MAX}" "${OVERLAP}" "${NSMEAR}" \
  --beam-energy-gev "${BEAM_ENERGY_GEV_ARG}" \
  "${SMEAR_GEOMETRY_ARGS[@]}" \
  --input-preview-pdf "${INPUT_PREVIEW_PDF}" \
  --confirm-inputs
SMEAR_STATUS=$?
set -e
if [[ ${SMEAR_STATUS} -eq 20 ]]; then
  echo "Pipeline stopped before optimization because the input histograms were not approved."
  exit 0
fi
if [[ ${SMEAR_STATUS} -ne 0 ]]; then
  echo "[ERROR] Smearing fitter exited with status ${SMEAR_STATUS}." >&2
  exit "${SMEAR_STATUS}"
fi

if [[ ! -f "${GENERATED_SECTION_MAP_FILE}" ]]; then
  echo "[ERROR] Section map was not generated: ${GENERATED_SECTION_MAP_FILE}" >&2
  exit 1
fi
if [[ ! -f "${GENERATED_INTERP_FILE}" ]]; then
  echo "[ERROR] Interpolated map ROOT was not generated: ${GENERATED_INTERP_FILE}" >&2
  exit 1
fi

if [[ "${SECTION_MAP_FILE}" != "${GENERATED_SECTION_MAP_FILE}" ]]; then
  cp -f "${GENERATED_SECTION_MAP_FILE}" "${SECTION_MAP_FILE}"
fi
if [[ "${INTERP_FILE}" != "${GENERATED_INTERP_FILE}" ]]; then
  cp -f "${GENERATED_INTERP_FILE}" "${INTERP_FILE}"
fi

run_simc_analysis "Step 4/4: Regenerating smeared simulation ROOT output" "section"

echo ""
echo "============================================================================"
echo "Pipeline complete"
echo "  kin:             ${KIN}"
echo "  target:          ${TARGET}"
echo "  mode hint:       ${MODE_HINT}"
echo "  smear rng seed:  ${SMEAR_RANDOM_SEED}"
echo "  beam energy GeV: ${BEAM_ENERGY_GEV_ARG}"
echo "  nps z [cm]:      ${NPS_Z_NPS_CM_ARG}"
if [[ -n "${NPS_THETA_NPS_DEG_ARG}" ]]; then
  echo "  nps theta [deg]: ${NPS_THETA_NPS_DEG_ARG} (NPS->Hall)"
else
  echo "  nps angle [deg]: ${NPS_ANGLE_DEG_ARG} (physical; fitter uses negative)"
fi
echo "  acceptance cfg:  ${ACCEPTANCE_CONFIG}"
echo "  combined input:  ${COMBINED_FILE}"
echo "  sim input:       ${SIM_FILE}"
echo "  smear fit root:  ${OUT_FILE}"
echo "  section map:     ${SECTION_MAP_FILE}"
echo "  interp map:      ${INTERP_FILE}"
echo "  input preview:   ${INPUT_PREVIEW_PDF}"
echo "  sim smeared out: ${SIM_SMEARED_FILE}"
echo "============================================================================"
