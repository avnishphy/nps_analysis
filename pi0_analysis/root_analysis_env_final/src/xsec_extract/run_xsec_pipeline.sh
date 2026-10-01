#!/usr/bin/env bash
# ============================================================================
# NPS pi0 xsec extraction pipeline wrapper
# - Resolves canonical input/output paths per kinematic setting
# - Selects and compiles either supported cross-section extractor
# - Runs extraction with validated paths
#
# Commands for x36_5_407 and x36_4:
# ./src/xsec_extract/run_xsec_pipeline.sh --kin KinC_x36_5_407 --target LH2 --root-dir output/KinC_x36_5_407/KinC_x36_5/root/ --sim-file /volatile/hallc/nps/singhav/nps_smearing/smear_x36_5_407/smearing_output/KinC_x36_5_407/root/simc_pi0_analysis_output_smeared.root --vertex_simc_file output/simc/simc_x36_5_407/worksim/ --mmiss-lower 0.6 --mmiss-upper 1.1 --partons --positive-xsec --mmiss_select ellipse --xsec_config xsec_config_x36_5_407.json --fit-objective scaled-poisson
# ./src/xsec_extract/run_xsec_pipeline.sh --kin KinC_x36_4 --target LH2 --root-dir output/KinC_x36_4/root/ --sim-file /volatile/hallc/nps/singhav/nps_smearing/smear_x36_4/smearing_output/KinC_x36_4/root/simc_pi0_analysis_output_smeared.root --vertex_simc_file output/simc/nps_simc_20260824_135058/worksim/simc_gfortran_updated/worksim/ --mmiss-lower 0.6 --mmiss-upper 1.1 --partons --positive-xsec --mmiss_select ellipse --xsec_config xsec_config_x36_4.json --fit-objective scaled-poisson

# src/xsec_extract/run_xsec_pipeline.sh \
#   --kin KinC_x36_4 \
#   --data-file output/KinC_x36_4/root/combined_branches_LH2.root \
#   --sim-file /volatile/hallc/nps/singhav/nps_smearing/smear_x36_4/smearing_output/KinC_x36_4/root/simc_pi0_analysis_output_smeared.root \
#   --vertex_simc_file output/simc/nps_simc_20260824_135058/worksim/simc_gfortran_updated/worksim/nps_excl_pi0_x36_4.root \
#   --mmiss_select ellipse --mmiss-cut-file output/KinC_x36_4/root/combined_branches_LH2_combined_2d_mass_cut_debug.txt --mmiss-lower 0.6 --mmiss-upper 1.1 --prepare-forward-inputs --out-dir output/KinC_x36_4/forward_cache'
# 
# Default extraction uses the missing-mass window; --mmiss_select chooses a
# stored combined-data selector. Full 0-2.5 GeV spectra remain diagnostic.
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

ROOT_CMD="${ROOT_CMD:-root}"
CXX_CMD="${CXX:-g++}"
XSEC_METHOD="${NPS_XSEC_METHOD:-no-simc-model}"
XSEC_SRC=""
XSEC_CONFIG="${NPS_XSEC_CONFIG:-}"

KIN="${NPS_XSEC_KIN:-}"
TARGET="${NPS_XSEC_TARGET:-LH2}"
OUTPUT_BASE="${NPS_OUTPUT_BASE:-${REPO_ROOT}/output}"
ROOT_DIR="${NPS_XSEC_ROOT_DIR:-}"

DATA_FILE="${NPS_XSEC_DATA_FILE:-}"
SIM_FILE="${NPS_XSEC_SIM_FILE:-}"
VERTEX_SIMC_FILE="${NPS_XSEC_VERTEX_SIMC_FILE:-}"

OUT_DIR="${NPS_XSEC_OUT_DIR:-}"
OUT_ROOT="${NPS_XSEC_OUT_ROOT:-}"
OUT_CSV="${NPS_XSEC_OUT_CSV:-}"
OUT_SLICE_CSV="${NPS_XSEC_OUT_SLICE_CSV:-}"
ALL_PLOTS_PDF="${NPS_XSEC_ALL_PLOTS_PDF:-}"
# Empty means use the selected extractor's configuration header. Explicit
# environment variables or CLI options below still take precedence.
MMISS_LOWER="${NPS_XSEC_MMISS_LOWER:-}"
MMISS_UPPER="${NPS_XSEC_MMISS_UPPER:-}"
MMISS_SELECT="${NPS_XSEC_MMISS_SELECT:-}"
MMISS_CUT_FILE="${NPS_XSEC_MMISS_CUT_FILE:-}"
TARGET_CONTAM="${NPS_XSEC_TARGET_CONTAM:-}"
TARGET_CONTAM_ERR="${NPS_XSEC_TARGET_CONTAM_ERR:-}"
NORMALIZE_MMISS="${NPS_XSEC_NORMALIZE_MMISS:-0}"
MODEL_ID=""
SIMC_YIELD_SCALE=""
SIMC_EBEAM=""
MODEL_FIXED=0
MODEL_FREE=""
MODEL_MAX_ITERATIONS=""
MODEL_MAX_EVALUATIONS=""
MODEL_TOLERANCE=""

# Numerical controls diagnose response rank and MC-variance convergence.
# They never rescale the response to match the measured data integral.
SVD_RANK_TOLERANCE=1e-10
MC_MAX_ITERATIONS=100
MC_FIT_TOLERANCE=1e-6
FIT_VARIANCE="${NPS_XSEC_FIT_VARIANCE:-finite-mc}"
FIT_OBJECTIVE="${NPS_XSEC_FIT_OBJECTIVE:-gaussian}"
SCALED_EMPTY_SCALE="${NPS_XSEC_SCALED_EMPTY_SCALE:-auto}"
POSITIVE_XSEC="${NPS_XSEC_POSITIVE_XSEC:-0}"
POSITIVE_XSEC_EXPLICIT=0
PARTONS_ENABLED=0
PARTONS_WARMUPS="${NPS_PARTONS_WARMUPS:-10000}"
PARTONS_CALLS="${NPS_PARTONS_CALLS:-100000}"
SOFTWARE_ROOT="${NPS_SOFTWARE_ROOT:-/group/nps/singhav/software}"
PARTONS_ROOT="${NPS_PARTONS_ROOT:-${SOFTWARE_ROOT}/partons}"

QUIET=0
NO_DIAGNOSTICS=0
NO_PDF=0
NO_PNG=0
PREPARE_FORWARD_INPUTS=0

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

print_help() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --prepare-forward-inputs    Export selected events before rectangular binning; no fit/plots
  --xsec-method <method>       no-simc-model (default) or simc-model
  --simc-model                 Shortcut for --xsec-method simc-model
  --no-simc-model              Shortcut for --xsec-method no-simc-model
  --xsec_config <json>        JSON preset in xsec_config/ (default: inferred from --kin)
  --kin <Kin_old>             Kinematic setting (recommended)
  --target <name>             Combined target token (default: LH2)
  --output-base <path>        Output base directory (default: repo/output)
  --root-dir <path>           Root directory with combined/sim files
  --data-file <path>          Combined data ROOT file
  --sim-file <path>           Simulation ROOT file
  --vertex_simc_file <file|dir> Required for no-simc-model: original exclusive h10 ROOT or worksim directory
  --out-dir <path>            Xsec output directory
  --out-root <path>           Xsec output ROOT file
  --out-csv <path>            Xsec output summary CSV
  --out-slice-csv <path>      Xsec output slice CSV
  --all-plots-pdf <path>      Xsec combined plots PDF
  --mmiss-lower <GeV>         Shared data/exclusive-SIMC lower bound (default: extractor config)
  --mmiss-upper <GeV>         Shared data/exclusive-SIMC upper bound (default: extractor config)
  --mmiss_select <mode>       mcd or ellipse (no-simc-model only)
  --mmiss-cut-file <path>     Combined-data geometric cut metadata (no-simc-model only)
  --target-contam <factor>    Data yield divisor (default: extractor config; use 1 to omit)
  --target-contam-err <factor> Absolute divisor uncertainty (default: extractor config)
  --normalize_mmiss          Historical selected-yield rescaling, permitted
                             only with --fixed-default-model; removes absolute sensitivity
  --normalize-simc-to-data  Clear-name alias for --normalize_mmiss
  --simc-yield-scale <x>     Independently justified SIMC normalization (default: 1)
  --ebeam <GeV>              Beam energy for reporting or simc-model epsilon fallback
  --model <id>               SIMC physics-model identifier
  --fixed-default-model      SIMC model: no parameter fit
  --model-free <names>       SIMC model: comma-separated Fortran coefficients
  --model-max-iterations <n> --model-max-evaluations <n>
  --model-tolerance <x>      SIMC model Minuit2 controls
  Bin edges and optional diamond: edit the selected xsec_config/xsec_config*.json file
  --svd-rank-tolerance <float> Relative singular-value rank cutoff (default: 1e-10)
  --mc-max-iterations <int>    MC-variance fit iteration limit (default: 30)
  --mc-fit-tolerance <float>   Relative MC-variance fit convergence (default: 1e-6)
  --fit-variance <data|finite-mc> Data-only Eq. 5.23 or finite-MC extension (default)
  --fit-objective <gaussian|scaled-poisson>  Fit statistic (default gaussian)
  --scaled-empty-scale <auto|slice|global>  Empty-row scale sensitivity
  --positive-xsec             Require full angular cross section >= 0; LT/TT signed
  --no-positive-xsec          Disable positivity (default; overrides environment)
  --partons                  Native GK06/GPDGK19 pi0 projection and comparison plots
  --partons-warmups <int>    DVMP CFF MC warm-up calls (default: 10000)
  --partons-calls <int>      DVMP CFF MC calls (default: 100000)
  --quiet                     Pass --quiet into xsec executable
  --no-diagnostics            Omit diagnostic plots; retain core fit plots
  --no-pdf                    Omit individual and combined PDFs
  --no-png                    Omit PNG plots
  --help                      Show this message

Environment overrides:
  ROOT_CMD, CXX, NPS_XSEC_METHOD, NPS_XSEC_CONFIG,
  NPS_XSEC_KIN, NPS_XSEC_TARGET, NPS_OUTPUT_BASE, NPS_XSEC_ROOT_DIR,
  NPS_XSEC_DATA_FILE, NPS_XSEC_SIM_FILE, NPS_XSEC_VERTEX_SIMC_FILE,
  NPS_XSEC_OUT_DIR, NPS_XSEC_OUT_ROOT, NPS_XSEC_OUT_CSV,
  NPS_XSEC_OUT_SLICE_CSV, NPS_XSEC_ALL_PLOTS_PDF.
  NPS_XSEC_MMISS_LOWER, NPS_XSEC_MMISS_UPPER,
  NPS_XSEC_MMISS_SELECT, NPS_XSEC_MMISS_CUT_FILE,
  NPS_XSEC_TARGET_CONTAM, NPS_XSEC_TARGET_CONTAM_ERR,
  NPS_XSEC_NORMALIZE_MMISS (0 or 1),
  NPS_XSEC_POSITIVE_XSEC (0 or 1; last explicit on/off flag wins),
  NPS_PARTONS_WARMUPS, NPS_PARTONS_CALLS, NPS_SOFTWARE_ROOT, NPS_PARTONS_ROOT.

Defaults when --kin is provided:
  root-dir      = <output-base>/<sanitize(kin)>/root
  data-file     = <root-dir>/combined_branches_<sanitize(target)>.root
  sim-file      = <root-dir>/simc_pi0_analysis_output_smeared.root
  no-simc-model: <output-base>/<kin>/xsec with *_no_simc_model* filenames
  simc-model:    <output-base>/<kin>/xsec_simc_model with *_simc_model* filenames
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --xsec-method)
      XSEC_METHOD="$2"
      shift 2
      ;;
    --xsec_config|--xsec-config)
      if [[ $# -lt 2 || -z "$2" || "$2" == --* ]]; then
        echo "[ERROR] --xsec_config requires a JSON preset name or path." >&2
        exit 1
      fi
      XSEC_CONFIG="$2"
      shift 2
      ;;
    --simc-model)
      XSEC_METHOD="simc-model"
      shift
      ;;
    --no-simc-model)
      XSEC_METHOD="no-simc-model"
      shift
      ;;
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
    --data-file)
      DATA_FILE="$2"
      shift 2
      ;;
    --sim-file)
      SIM_FILE="$2"
      shift 2
      ;;
    --vertex_simc_file|--vertex-simc-file)
      VERTEX_SIMC_FILE="$2"
      shift 2
      ;;
    --out-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --out-root)
      OUT_ROOT="$2"
      shift 2
      ;;
    --out-csv)
      OUT_CSV="$2"
      shift 2
      ;;
    --out-slice-csv)
      OUT_SLICE_CSV="$2"
      shift 2
      ;;
    --all-plots-pdf)
      ALL_PLOTS_PDF="$2"
      shift 2
      ;;
    --mmiss-lower)
      MMISS_LOWER="$2"
      shift 2
      ;;
    --mmiss-upper)
      MMISS_UPPER="$2"
      shift 2
      ;;
    --mmiss_select|--mmiss-select)
      MMISS_SELECT="$2"
      shift 2
      ;;
    --mmiss-cut-file)
      MMISS_CUT_FILE="$2"
      shift 2
      ;;
    --target-contam)
      TARGET_CONTAM="$2"
      shift 2
      ;;
    --target-contam-err)
      TARGET_CONTAM_ERR="$2"
      shift 2
      ;;
    --normalize_mmiss|--normalize-mmiss|--normalize-simc-to-data)
      NORMALIZE_MMISS=1
      shift
      ;;
    --simc-yield-scale) SIMC_YIELD_SCALE="$2"; shift 2 ;;
    --ebeam) SIMC_EBEAM="$2"; shift 2 ;;
    --model) MODEL_ID="$2"; shift 2 ;;
    --fixed-default-model) MODEL_FIXED=1; shift ;;
    --model-free) MODEL_FREE="$2"; shift 2 ;;
    --model-max-iterations) MODEL_MAX_ITERATIONS="$2"; shift 2 ;;
    --model-max-evaluations) MODEL_MAX_EVALUATIONS="$2"; shift 2 ;;
    --model-tolerance) MODEL_TOLERANCE="$2"; shift 2 ;;

    --svd-rank-tolerance)
      SVD_RANK_TOLERANCE="$2"
      shift 2
      ;;
    --mc-max-iterations)
      MC_MAX_ITERATIONS="$2"
      shift 2
      ;;
    --mc-fit-tolerance)
      MC_FIT_TOLERANCE="$2"
      shift 2
      ;;
    --fit-variance)
      FIT_VARIANCE="$2"
      shift 2
      ;;
    --fit-objective)
      FIT_OBJECTIVE="$2"
      shift 2
      ;;
    --prepare-forward-inputs)
      PREPARE_FORWARD_INPUTS=1
      NO_DIAGNOSTICS=1
      NO_PDF=1
      NO_PNG=1
      shift
      ;;
    --scaled-empty-scale)
      SCALED_EMPTY_SCALE="$2"
      shift 2
      ;;
    --positive-xsec)
      POSITIVE_XSEC_EXPLICIT=1
      POSITIVE_XSEC=1
      shift
      ;;
    --no-positive-xsec)
      POSITIVE_XSEC_EXPLICIT=1
      POSITIVE_XSEC=0
      shift
      ;;
    --partons)
      PARTONS_ENABLED=1
      shift
      ;;
    --partons-warmups)
      PARTONS_WARMUPS="$2"
      shift 2
      ;;
    --partons-calls)
      PARTONS_CALLS="$2"
      shift 2
      ;;
    --quiet)
      QUIET=1
      shift
      ;;
    --no-diagnostics)
      NO_DIAGNOSTICS=1
      shift
      ;;
    --no-pdf)
      NO_PDF=1
      shift
      ;;
    --no-png)
      NO_PNG=1
      shift
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

case "${XSEC_METHOD}" in
  no-simc-model|no_simc_model)
    XSEC_METHOD="no-simc-model"
    XSEC_SRC="${SCRIPT_DIR}/excl_xsec_pi0_analysis_no_simc_model.C"
    XSEC_BASENAME="excl_xsec_pi0_analysis_no_simc_model"
    XSEC_OUT_SUBDIR="xsec"
    XSEC_FALLBACK_DIR="output_pi0_xsec_no_simc_model"
    ;;
  simc-model|simc_model)
    XSEC_METHOD="simc-model"
    XSEC_SRC="${SCRIPT_DIR}/excl_xsec_pi0_analysis_simc_model.C"
    XSEC_BASENAME="excl_xsec_pi0_analysis_simc_model"
    XSEC_OUT_SUBDIR="xsec_simc_model"
    XSEC_FALLBACK_DIR="output_pi0_xsec_simc_model"
    ;;
  *)
    echo "[ERROR] Invalid --xsec-method '${XSEC_METHOD}'; use no-simc-model or simc-model." >&2
    exit 1
    ;;
esac

CONFIG_DIR="$(cd "${SCRIPT_DIR}/xsec_config" && pwd -P)"
if [[ -z "${XSEC_CONFIG}" && -n "${KIN}" ]]; then
  candidate="xsec_config_${KIN#KinC_}.json"
  if [[ -f "${CONFIG_DIR}/${candidate}" ]]; then XSEC_CONFIG="${candidate}"; fi
fi
if [[ -z "${XSEC_CONFIG}" ]]; then
  echo "[ERROR] Select a JSON preset in ${CONFIG_DIR} with --xsec_config." >&2
  exit 1
fi
case "${XSEC_CONFIG}" in
  /*) CONFIG_CANDIDATE="${XSEC_CONFIG}" ;;
  xsec_config/*) CONFIG_CANDIDATE="${SCRIPT_DIR}/${XSEC_CONFIG}" ;;
  */*) CONFIG_CANDIDATE="${REPO_ROOT}/${XSEC_CONFIG}" ;;
  *) CONFIG_CANDIDATE="${CONFIG_DIR}/${XSEC_CONFIG}" ;;
esac
if [[ ! -f "${CONFIG_CANDIDATE}" ]]; then
  echo "[ERROR] Xsec JSON preset not found: ${CONFIG_CANDIDATE}" >&2
  exit 1
fi
XSEC_CONFIG_FILE="$(readlink -f "${CONFIG_CANDIDATE}")"
CONFIG_NAME="$(basename "${XSEC_CONFIG_FILE}")"
if [[ "$(dirname "${XSEC_CONFIG_FILE}")" != "${CONFIG_DIR}" ||
      ! "${CONFIG_NAME}" =~ ^xsec_config[A-Za-z0-9_]*\.json$ ]]; then
  echo "[ERROR] --xsec_config must select an xsec_config*.json file in ${CONFIG_DIR}." >&2
  exit 1
fi

if [[ "${XSEC_METHOD}" == "simc-model" && ( -n "${MMISS_SELECT}" || -n "${MMISS_CUT_FILE}" ) ]]; then
  echo "[ERROR] --mmiss_select and --mmiss-cut-file apply only to --xsec-method no-simc-model." >&2
  exit 1
fi
case "${MMISS_SELECT}" in
  ""|mcd|ellipse) ;;
  *) echo "[ERROR] --mmiss_select must be mcd or ellipse." >&2; exit 1 ;;
esac
if [[ "${XSEC_METHOD}" == "simc-model" && "${POSITIVE_XSEC}" == "1" ]]; then
  echo "[ERROR] --positive-xsec is supported only with --xsec-method no-simc-model." >&2
  exit 1
fi
if [[ "${FIT_OBJECTIVE}" != "gaussian" && "${FIT_OBJECTIVE}" != "scaled-poisson" ]]; then
  echo "[ERROR] --fit-objective must be gaussian or scaled-poisson." >&2
  exit 1
fi
if [[ "${FIT_OBJECTIVE}" == "scaled-poisson" ]]; then
  if [[ "${POSITIVE_XSEC_EXPLICIT}" == "1" && "${POSITIVE_XSEC}" == "0" ]]; then
    echo "[ERROR] scaled-poisson requires angular positivity; remove --no-positive-xsec." >&2
    exit 1
  fi
  POSITIVE_XSEC=1
fi
if [[ "${XSEC_METHOD}" == "simc-model" && "${FIT_OBJECTIVE}" != "gaussian" ]]; then
  echo "[ERROR] --fit-objective scaled-poisson applies only to no-simc-model." >&2
  exit 1
fi
if [[ "${FIT_VARIANCE}" != "data" && "${FIT_VARIANCE}" != "finite-mc" ]]; then
  echo "[ERROR] --fit-variance must be data or finite-mc." >&2
  exit 1
fi
if [[ "${XSEC_METHOD}" == "simc-model" && "${FIT_VARIANCE}" != "finite-mc" ]]; then
  echo "[ERROR] --fit-variance applies only to --xsec-method no-simc-model." >&2
  exit 1
fi

case "${NORMALIZE_MMISS}" in
  0|1) ;;
  *)
    echo "[ERROR] NPS_XSEC_NORMALIZE_MMISS must be 0 or 1." >&2
    exit 1
    ;;
esac
if [[ "${NORMALIZE_MMISS}" == "1" && "${XSEC_METHOD}" == "simc-model" && "${MODEL_FIXED}" != "1" ]]; then
  echo "[ERROR] --normalize_mmiss conflicts with iterative absolute model fitting; use --fixed-default-model." >&2
  exit 1
fi
if [[ "${NORMALIZE_MMISS}" == "1" ]]; then
  if [[ "${XSEC_METHOD}" != "simc-model" ]]; then
    echo "[ERROR] --normalize_mmiss is currently supported only with --xsec-method simc-model." >&2
    exit 1
  fi
  # Historical flag semantics: the data integral sets the SIMC scale and
  # replaces, rather than compounds, the separate target correction.
  TARGET_CONTAM=1.0
  TARGET_CONTAM_ERR=0.0
fi

case "${POSITIVE_XSEC}" in
  0|1) ;;
  *)
    echo "[ERROR] NPS_XSEC_POSITIVE_XSEC must be 0 or 1 (or override with --positive-xsec/--no-positive-xsec)." >&2
    exit 1
    ;;
esac

KIN="$(trim_ws "${KIN}")"
TARGET="$(trim_ws "${TARGET}")"
OUTPUT_BASE="$(to_abs_path "${OUTPUT_BASE}")"

if [[ -n "${KIN}" && -z "${ROOT_DIR}" ]]; then
  ROOT_DIR="${OUTPUT_BASE}/$(sanitize_name "${KIN}")/root"
fi
ROOT_DIR="$(to_abs_path "${ROOT_DIR}")"

if [[ -n "${ROOT_DIR}" && -z "${DATA_FILE}" ]]; then
  DATA_FILE="${ROOT_DIR}/combined_branches_$(sanitize_name "${TARGET}").root"
fi
if [[ -n "${ROOT_DIR}" && -z "${SIM_FILE}" ]]; then
  SIM_FILE="${ROOT_DIR}/simc_pi0_analysis_output_smeared.root"
fi

if [[ -n "${KIN}" && -z "${OUT_DIR}" ]]; then
  OUT_DIR="${OUTPUT_BASE}/$(sanitize_name "${KIN}")/${XSEC_OUT_SUBDIR}"
elif [[ -n "${ROOT_DIR}" && -z "${OUT_DIR}" ]]; then
  OUT_DIR="$(dirname "${ROOT_DIR}")/${XSEC_OUT_SUBDIR}"
fi
if [[ -z "${OUT_DIR}" ]]; then
  OUT_DIR="${REPO_ROOT}/${XSEC_FALLBACK_DIR}"
fi
OUT_DIR="$(to_abs_path "${OUT_DIR}")"

if [[ -z "${OUT_ROOT}" ]]; then
  OUT_ROOT="${OUT_DIR}/${XSEC_BASENAME}_output.root"
fi
if [[ -z "${OUT_CSV}" ]]; then
  OUT_CSV="${OUT_DIR}/${XSEC_BASENAME}_summary.csv"
fi
if [[ -z "${OUT_SLICE_CSV}" ]]; then
  OUT_SLICE_CSV="${OUT_DIR}/${XSEC_BASENAME}_slice_summary.csv"
fi
if [[ -z "${ALL_PLOTS_PDF}" ]]; then
  ALL_PLOTS_PDF="${OUT_DIR}/all_generated_plots_${XSEC_METHOD//-/_}.pdf"
fi

DATA_FILE="$(to_abs_path "${DATA_FILE}")"
SIM_FILE="$(to_abs_path "${SIM_FILE}")"
VERTEX_SIMC_FILE="$(to_abs_path "${VERTEX_SIMC_FILE}")"
OUT_ROOT="$(to_abs_path "${OUT_ROOT}")"
OUT_CSV="$(to_abs_path "${OUT_CSV}")"
OUT_SLICE_CSV="$(to_abs_path "${OUT_SLICE_CSV}")"
ALL_PLOTS_PDF="$(to_abs_path "${ALL_PLOTS_PDF}")"
MMISS_CUT_FILE="$(to_abs_path "${MMISS_CUT_FILE}")"

if [[ -z "${DATA_FILE}" || -z "${SIM_FILE}" ]]; then
  echo "[ERROR] Unable to resolve data/sim input files." >&2
  echo "        Provide --data-file and --sim-file, or use --kin/--root-dir." >&2
  exit 1
fi
if [[ ! -f "${DATA_FILE}" ]]; then
  echo "[ERROR] Data file not found: ${DATA_FILE}" >&2
  exit 1
fi
if [[ ! -f "${SIM_FILE}" ]]; then
  echo "[ERROR] Simulation file not found: ${SIM_FILE}" >&2
  exit 1
fi
if [[ -n "${VERTEX_SIMC_FILE}" && ! -f "${VERTEX_SIMC_FILE}" && ! -d "${VERTEX_SIMC_FILE}" ]]; then
  echo "[ERROR] Original SIMC file/directory not found: ${VERTEX_SIMC_FILE}" >&2
  exit 1
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
if [[ ! -f "${XSEC_SRC}" ]]; then
  echo "[ERROR] Missing xsec source file: ${XSEC_SRC}" >&2
  exit 1
fi

if [[ "${PREPARE_FORWARD_INPUTS}" -eq 1 ]]; then
  if [[ "${XSEC_METHOD}" != "no-simc-model" || "${PARTONS_ENABLED}" -eq 1 ]]; then
    echo "[ERROR] Forward event preparation requires no-simc-model without PARTONS." >&2
    exit 1
  fi
  for cache_name in data_events.csv mc_events.csv forward_cache_manifest.json; do
    if [[ -e "${OUT_DIR}/${cache_name}" || -e "${OUT_DIR}/${cache_name}.partial" ]]; then
      echo "[ERROR] Forward event cache already exists or is partial: ${OUT_DIR}/${cache_name}" >&2
      exit 1
    fi
  done
fi
mkdir -p "${OUT_DIR}" "$(dirname "${OUT_ROOT}")" "$(dirname "${OUT_CSV}")" "$(dirname "${OUT_SLICE_CSV}")" "$(dirname "${ALL_PLOTS_PDF}")"
if [[ "${MMISS_SELECT}" == "mcd" || "${MMISS_SELECT}" == "ellipse" ]]; then
  if [[ -z "${MMISS_CUT_FILE}" ]]; then
    candidate="${DATA_FILE%.*}_combined_2d_mass_cut_debug.txt"
    key="cov_mpi0_mpi0"
    if [[ "${MMISS_SELECT}" == "mcd" ]]; then key="mcd_cov_mpi0_mpi0"; fi
    if [[ -f "${candidate}" ]] && grep -q "^${key}=" "${candidate}" && grep -q '^mpi0_min=' "${candidate}"; then
      MMISS_CUT_FILE="${candidate}"
    else
      MMISS_CUT_FILE="${OUT_DIR}/$(basename "${DATA_FILE%.*}")_combined_2d_mass_cut_debug.txt"
      echo "[mass-cut] Exporting verified geometry from ${DATA_FILE}"
      python3 "${SCRIPT_DIR}/../analysis/export_combined_mass_cut_metadata.py" "${DATA_FILE}" --out "${MMISS_CUT_FILE}"
    fi
  fi
fi

echo "============================================================================"
echo "NPS pi0 xsec pipeline"
echo "  xsec method:     ${XSEC_METHOD}"
echo "  xsec config:     ${XSEC_CONFIG_FILE}"
echo "  kin:             ${KIN:-<not-set>}"
echo "  target:          ${TARGET}"
echo "  output base:     ${OUTPUT_BASE}"
echo "  root dir:        ${ROOT_DIR:-<not-set>}"
echo "  data file:       ${DATA_FILE}"
echo "  sim file:        ${SIM_FILE}"
if [[ "${XSEC_METHOD}" == "no-simc-model" ]]; then
  echo "  vertex SIMC:     ${VERTEX_SIMC_FILE:-<off>}"
else
  echo "  SIMC weighting:  simc_yield_scale*(full_weight/sigcm)*model(vertex,p)"
fi
echo "  selection:       ${MMISS_SELECT:-window} (${MMISS_LOWER:-<config>} < Mmiss < ${MMISS_UPPER:-<config>} GeV for window)"
echo "  target factor:   ${TARGET_CONTAM:-<config>} +/- ${TARGET_CONTAM_ERR:-<config>} (data divided by factor)"
echo "  SIMC/data yield normalization: ${NORMALIZE_MMISS} (global selected-yield area match)"
if [[ "${XSEC_METHOD}" == "no-simc-model" ]]; then
  echo "  SVD rank tol:    ${SVD_RANK_TOLERANCE}"
  echo "  MC iterations:   ${MC_MAX_ITERATIONS} (convergence ${MC_FIT_TOLERANCE})"
  echo "  angular positivity: ${POSITIVE_XSEC} (zero allowed; LT/TT remain signed)"
fi
echo "  PARTONS:         ${PARTONS_ENABLED} (warmups ${PARTONS_WARMUPS}, calls ${PARTONS_CALLS})"
echo "  plots:           PDF=$((1-NO_PDF)) PNG=$((1-NO_PNG)) diagnostics=$((1-NO_DIAGNOSTICS))"
echo "  out dir:         ${OUT_DIR}"
echo "  out root:        ${OUT_ROOT}"
echo "  out csv:         ${OUT_CSV}"
echo "  out slice csv:   ${OUT_SLICE_CSV}"
echo "  all plots pdf:   ${ALL_PLOTS_PDF}"
echo "============================================================================"

BUILD_DIR="$(mktemp -d "${OUT_DIR}/.build_xsec.XXXXXX")"
trap 'rm -rf "${BUILD_DIR}"' EXIT
XSEC_BIN="${BUILD_DIR}/${XSEC_BASENAME}"

echo "[build] Generating xsec_config.h from ${XSEC_CONFIG_FILE}"
python3 "${SCRIPT_DIR}/generate_xsec_config.py" "${XSEC_CONFIG_FILE}" "${BUILD_DIR}/xsec_config.h"
echo "[build] Compiling xsec executable"
if [[ "${PARTONS_ENABLED}" -eq 1 ]]; then
  # The installed v5 C++ example supplies a verified library/include layout.
  # Keep PARTONS optional so xsec jobs on hosts without this installation
  # retain the original ROOT-only build. Do not silently substitute PDF sets.
  PARTONS_HEADER="${SCRIPT_DIR}/partons_pi0_projection.h"
  PARTONS_LIBS=(
    "${PARTONS_ROOT}/lib64/libsfml-system.so"
    "${PARTONS_ROOT}/lib/libcln.so"
    "${PARTONS_ROOT}/lib/libElementaryUtils.so"
    "${PARTONS_ROOT}/lib/libNumA++.so"
    "${PARTONS_ROOT}/lib/libPARTONS.so"
    "${SOFTWARE_ROOT}/python/lib/libgsl.so"
    "${SOFTWARE_ROOT}/python/lib/libgslcblas.so"
    "${SOFTWARE_ROOT}/apfel/lib/libapfelxx.so"
    "${SOFTWARE_ROOT}/lhapdf/lib/libLHAPDF.so"
    "/usr/lib64/libxml2.so"
    # The installed LHAPDF was built against CXXABI_1.3.15. Explicitly
    # resolve that symbol from its bundled C++ runtime at link time; the
    # same directory is first in LD_LIBRARY_PATH during model evaluation.
    "${SOFTWARE_ROOT}/lhapdf/lib/libstdc++.so"
  )
  [[ -f "${PARTONS_HEADER}" ]] || { echo "[ERROR] Missing PARTONS adapter: ${PARTONS_HEADER}" >&2; exit 1; }
  PARTONS_SCHEMA="${SOFTWARE_ROOT}/src/partons/partons-example/data/xmlSchema.xsd"
  [[ -f "${PARTONS_SCHEMA}" ]] || { echo "[ERROR] Missing PARTONS XML schema: ${PARTONS_SCHEMA}" >&2; exit 1; }
  for library in "${PARTONS_LIBS[@]}"; do
    [[ -f "${library}" ]] || { echo "[ERROR] Missing PARTONS dependency: ${library}" >&2; exit 1; }
  done
  # PARTONS::init locates partons.properties beside argv[0]. Both properties
  # files are generated next to this temporary executable and removed with
  # it after the run. Absolute logger/schema paths work from any repo cwd.
  cat > "${BUILD_DIR}/partons.properties" <<EOF
log.file.path = ${BUILD_DIR}/logger.properties
xml.schema.file.path = ${PARTONS_SCHEMA}
computation.nb.processor = 1
gpd.service.batch.size = 1000
collinear_distribution.service.batch.size = 1000
ccf.service.batch.size = 1000
observable.service.batch.size = 1000
EOF
  cat > "${BUILD_DIR}/logger.properties" <<EOF
enable = true
default.level = WARN
print.mode = COUT
log.folder.path = ${BUILD_DIR}
EOF
  PARTONS_RPATH="${PARTONS_ROOT}/lib64:${PARTONS_ROOT}/lib:${SOFTWARE_ROOT}/python/lib:${SOFTWARE_ROOT}/apfel/lib:${SOFTWARE_ROOT}/lhapdf/lib"
  "${CXX_CMD}" "${XSEC_SRC}" -I"${BUILD_DIR}" -I"${SCRIPT_DIR}" $(root-config --cflags --libs) -lMinuit2 -O2 -std=c++17 \
    -DNPS_ENABLE_PARTONS -I"${PARTONS_ROOT}/include" \
    -I"${SOFTWARE_ROOT}/python/include" -I"${SOFTWARE_ROOT}/apfel/include" \
    -I"${SOFTWARE_ROOT}/lhapdf/include" -I/usr/include/libxml2 \
    -Wl,-rpath,"${PARTONS_RPATH}" "${PARTONS_LIBS[@]}" -o "${XSEC_BIN}"
  # The Conda LHAPDF install carries the libstdc++ runtime used by PARTONS.
  export LD_LIBRARY_PATH="${SOFTWARE_ROOT}/lhapdf/lib:${PARTONS_RPATH}:${LD_LIBRARY_PATH:-}"
else
  if [[ "${XSEC_METHOD}" == "simc-model" ]]; then
    "${CXX_CMD}" "${XSEC_SRC}" -I"${BUILD_DIR}" -I"${SCRIPT_DIR}" $(root-config --cflags --libs) -lMinuit2 -O2 -std=c++17 -o "${XSEC_BIN}"
  else
    "${CXX_CMD}" "${XSEC_SRC}" -I"${BUILD_DIR}" -I"${SCRIPT_DIR}" $(root-config --cflags --libs) -lMinuit2 -O2 -std=c++17 -o "${XSEC_BIN}"
  fi
fi
declare -a xsec_cmd=(
  "${XSEC_BIN}"
  --data-file "${DATA_FILE}"
  --sim-file "${SIM_FILE}"
  --out-dir "${OUT_DIR}"
  --out-root "${OUT_ROOT}"
  --out-csv "${OUT_CSV}"
  --out-slice-csv "${OUT_SLICE_CSV}"
  --all-plots-pdf "${ALL_PLOTS_PDF}"
)

if [[ -n "${MMISS_LOWER}" ]]; then xsec_cmd+=(--mmiss-lower "${MMISS_LOWER}"); fi
if [[ -n "${MMISS_UPPER}" ]]; then xsec_cmd+=(--mmiss-upper "${MMISS_UPPER}"); fi
if [[ -n "${TARGET_CONTAM}" ]]; then xsec_cmd+=(--target-contam "${TARGET_CONTAM}"); fi
if [[ -n "${TARGET_CONTAM_ERR}" ]]; then xsec_cmd+=(--target-contam-err "${TARGET_CONTAM_ERR}"); fi

if [[ "${XSEC_METHOD}" == "no-simc-model" ]]; then
  if [[ -n "${MMISS_SELECT}" ]]; then xsec_cmd+=(--mmiss_select "${MMISS_SELECT}"); fi
  if [[ "${PREPARE_FORWARD_INPUTS}" -eq 1 ]]; then xsec_cmd+=(--prepare-forward-inputs); fi
  if [[ -n "${MMISS_CUT_FILE}" ]]; then xsec_cmd+=(--mmiss-cut-file "${MMISS_CUT_FILE}"); fi
  xsec_cmd+=(
    --svd-rank-tolerance "${SVD_RANK_TOLERANCE}"
    --mc-max-iterations "${MC_MAX_ITERATIONS}"
    --mc-fit-tolerance "${MC_FIT_TOLERANCE}"
    --fit-variance "${FIT_VARIANCE}"
    --fit-objective "${FIT_OBJECTIVE}"
    --scaled-empty-scale "${SCALED_EMPTY_SCALE}"
  )
  if [[ "${POSITIVE_XSEC}" -eq 1 ]]; then
    xsec_cmd+=(--positive-xsec)
  else
    # Always forward the resolved choice so an environment value cannot
    # override an explicit --no-positive-xsec in the compiled executable.
    xsec_cmd+=(--no-positive-xsec)
  fi
  if [[ -n "${VERTEX_SIMC_FILE}" ]]; then
    xsec_cmd+=(--vertex_simc_file "${VERTEX_SIMC_FILE}")
  fi
elif [[ -n "${VERTEX_SIMC_FILE}" ]]; then
  echo "[WARN] vertex SIMC input ignored by simc-model extractor." >&2
fi

if [[ "${XSEC_METHOD}" == "simc-model" ]]; then
  if [[ -n "${SIMC_YIELD_SCALE}" ]]; then xsec_cmd+=(--simc-yield-scale "${SIMC_YIELD_SCALE}"); fi
  if [[ -n "${SIMC_EBEAM}" ]]; then xsec_cmd+=(--ebeam "${SIMC_EBEAM}"); fi
  if [[ -n "${MODEL_ID}" ]]; then xsec_cmd+=(--model "${MODEL_ID}"); fi
  if [[ "${MODEL_FIXED}" == "1" ]]; then xsec_cmd+=(--fixed-default-model); fi
  if [[ -n "${MODEL_FREE}" ]]; then xsec_cmd+=(--model-free "${MODEL_FREE}"); fi
  if [[ -n "${MODEL_MAX_ITERATIONS}" ]]; then xsec_cmd+=(--model-max-iterations "${MODEL_MAX_ITERATIONS}"); fi
  if [[ -n "${MODEL_MAX_EVALUATIONS}" ]]; then xsec_cmd+=(--model-max-evaluations "${MODEL_MAX_EVALUATIONS}"); fi
  if [[ -n "${MODEL_TOLERANCE}" ]]; then xsec_cmd+=(--model-tolerance "${MODEL_TOLERANCE}"); fi
fi
if [[ "${PARTONS_ENABLED}" -eq 1 ]]; then
  xsec_cmd+=(--partons --partons-warmups "${PARTONS_WARMUPS}" --partons-calls "${PARTONS_CALLS}")
fi
if [[ "${NORMALIZE_MMISS}" -eq 1 ]]; then
  xsec_cmd+=(--normalize_mmiss)
fi

if [[ -n "${KIN}" ]]; then
  xsec_cmd+=(--kin "${KIN}")
fi
if [[ -n "${ROOT_DIR}" ]]; then
  xsec_cmd+=(--root-dir "${ROOT_DIR}")
fi
xsec_cmd+=(--target "${TARGET}" --output-base "${OUTPUT_BASE}")

if [[ "${QUIET}" -eq 1 ]]; then
  xsec_cmd+=(--quiet)
fi
if [[ "${NO_DIAGNOSTICS}" -eq 1 ]]; then
  xsec_cmd+=(--no-diagnostics)
fi
if [[ "${NO_PDF}" -eq 1 ]]; then
  xsec_cmd+=(--no-pdf)
fi
if [[ "${NO_PNG}" -eq 1 ]]; then
  xsec_cmd+=(--no-png)
fi

echo "[run] Running xsec extraction"
if [[ "${PREPARE_FORWARD_INPUTS}" -eq 1 ]]; then
  cp "${BUILD_DIR}/xsec_config.h" "${OUT_DIR}/forward_export_config.h"
  python3 - "${OUT_DIR}" "${SCRIPT_DIR}" "${BUILD_DIR}/xsec_config.h" "${xsec_cmd[@]}" <<'PY'
import hashlib, json, os, pathlib, sys
out, source, config = map(pathlib.Path, sys.argv[1:4])
files = sorted(source.glob('*.h')) + [source/'xsec_config_template.h.in',
    source/'excl_xsec_pi0_analysis_no_simc_model.C', source/'run_xsec_pipeline.sh',
    source/'generate_xsec_config.py', config]
manifest = {'argv': sys.argv[4:], 'sources': [
    {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]}
temporary = out/'forward_export_provenance.json.partial'
with temporary.open('x') as stream:
    json.dump(manifest, stream, indent=2)
    stream.write('\n')
os.replace(temporary, out/'forward_export_provenance.json')
PY
fi
"${xsec_cmd[@]}"

echo "============================================================================"
echo "Xsec pipeline complete"
echo "============================================================================"
