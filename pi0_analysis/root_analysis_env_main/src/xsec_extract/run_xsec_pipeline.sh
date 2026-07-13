#!/usr/bin/env bash
# ============================================================================
# NPS pi0 xsec extraction pipeline wrapper
# - Resolves canonical input/output paths per kinematic setting
# - Compiles excl_xsec_pi0_analysis_no_simc_model.C
# - Runs extraction with validated paths
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

ROOT_CMD="${ROOT_CMD:-root}"
CXX_CMD="${CXX:-g++}"
XSEC_SRC="${SCRIPT_DIR}/excl_xsec_pi0_analysis_no_simc_model.C"

KIN="${NPS_XSEC_KIN:-}"
TARGET="${NPS_XSEC_TARGET:-LH2}"
OUTPUT_BASE="${NPS_OUTPUT_BASE:-${REPO_ROOT}/output}"
ROOT_DIR="${NPS_XSEC_ROOT_DIR:-}"

DATA_FILE="${NPS_XSEC_DATA_FILE:-}"
SIM_FILE="${NPS_XSEC_SIM_FILE:-}"

OUT_DIR="${NPS_XSEC_OUT_DIR:-}"
OUT_ROOT="${NPS_XSEC_OUT_ROOT:-}"
OUT_CSV="${NPS_XSEC_OUT_CSV:-}"
OUT_SLICE_CSV="${NPS_XSEC_OUT_SLICE_CSV:-}"
ALL_PLOTS_PDF="${NPS_XSEC_ALL_PLOTS_PDF:-}"

QUIET=0

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
  --kin <Kin_old>             Kinematic setting (recommended)
  --target <name>             Combined target token (default: LH2)
  --output-base <path>        Output base directory (default: repo/output)
  --root-dir <path>           Root directory with combined/sim files
  --data-file <path>          Combined data ROOT file
  --sim-file <path>           Simulation ROOT file
  --out-dir <path>            Xsec output directory
  --out-root <path>           Xsec output ROOT file
  --out-csv <path>            Xsec output summary CSV
  --out-slice-csv <path>      Xsec output slice CSV
  --all-plots-pdf <path>      Xsec combined plots PDF
  --quiet                     Pass --quiet into xsec executable
  --help                      Show this message

Environment overrides:
  ROOT_CMD, CXX,
  NPS_XSEC_KIN, NPS_XSEC_TARGET, NPS_OUTPUT_BASE, NPS_XSEC_ROOT_DIR,
  NPS_XSEC_DATA_FILE, NPS_XSEC_SIM_FILE,
  NPS_XSEC_OUT_DIR, NPS_XSEC_OUT_ROOT, NPS_XSEC_OUT_CSV,
  NPS_XSEC_OUT_SLICE_CSV, NPS_XSEC_ALL_PLOTS_PDF.

Defaults when --kin is provided:
  root-dir      = <output-base>/<sanitize(kin)>/root
  data-file     = <root-dir>/combined_branches_<sanitize(target)>.root
  sim-file      = <root-dir>/simc_pi0_analysis_output_smeared.root (fallback to unsmeared)
  out-dir       = <output-base>/<sanitize(kin)>/xsec
  out-root      = <out-dir>/excl_xsec_pi0_analysis_no_simc_model_output.root
  out-csv       = <out-dir>/excl_xsec_pi0_analysis_no_simc_model_summary.csv
  out-slice-csv = <out-dir>/excl_xsec_pi0_analysis_no_simc_model_slice_summary.csv
  all-plots-pdf = <out-dir>/all_generated_plots_no_simc_model.pdf
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
    --data-file)
      DATA_FILE="$2"
      shift 2
      ;;
    --sim-file)
      SIM_FILE="$2"
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
    --quiet)
      QUIET=1
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
  if [[ -f "${ROOT_DIR}/simc_pi0_analysis_output_smeared.root" ]]; then
    SIM_FILE="${ROOT_DIR}/simc_pi0_analysis_output_smeared.root"
  else
    SIM_FILE="${ROOT_DIR}/simc_pi0_analysis_output.root"
  fi
fi

if [[ -n "${KIN}" && -z "${OUT_DIR}" ]]; then
  OUT_DIR="${OUTPUT_BASE}/$(sanitize_name "${KIN}")/xsec"
elif [[ -n "${ROOT_DIR}" && -z "${OUT_DIR}" ]]; then
  OUT_DIR="$(dirname "${ROOT_DIR}")/xsec"
fi
if [[ -z "${OUT_DIR}" ]]; then
  OUT_DIR="${REPO_ROOT}/output_pi0_xsec_no_simc_model"
fi
OUT_DIR="$(to_abs_path "${OUT_DIR}")"

if [[ -z "${OUT_ROOT}" ]]; then
  OUT_ROOT="${OUT_DIR}/excl_xsec_pi0_analysis_no_simc_model_output.root"
fi
if [[ -z "${OUT_CSV}" ]]; then
  OUT_CSV="${OUT_DIR}/excl_xsec_pi0_analysis_no_simc_model_summary.csv"
fi
if [[ -z "${OUT_SLICE_CSV}" ]]; then
  OUT_SLICE_CSV="${OUT_DIR}/excl_xsec_pi0_analysis_no_simc_model_slice_summary.csv"
fi
if [[ -z "${ALL_PLOTS_PDF}" ]]; then
  ALL_PLOTS_PDF="${OUT_DIR}/all_generated_plots_no_simc_model.pdf"
fi

DATA_FILE="$(to_abs_path "${DATA_FILE}")"
SIM_FILE="$(to_abs_path "${SIM_FILE}")"
OUT_ROOT="$(to_abs_path "${OUT_ROOT}")"
OUT_CSV="$(to_abs_path "${OUT_CSV}")"
OUT_SLICE_CSV="$(to_abs_path "${OUT_SLICE_CSV}")"
ALL_PLOTS_PDF="$(to_abs_path "${ALL_PLOTS_PDF}")"

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

mkdir -p "${OUT_DIR}" "$(dirname "${OUT_ROOT}")" "$(dirname "${OUT_CSV}")" "$(dirname "${OUT_SLICE_CSV}")" "$(dirname "${ALL_PLOTS_PDF}")"

echo "============================================================================"
echo "NPS pi0 xsec pipeline"
echo "  kin:             ${KIN:-<not-set>}"
echo "  target:          ${TARGET}"
echo "  output base:     ${OUTPUT_BASE}"
echo "  root dir:        ${ROOT_DIR:-<not-set>}"
echo "  data file:       ${DATA_FILE}"
echo "  sim file:        ${SIM_FILE}"
echo "  out dir:         ${OUT_DIR}"
echo "  out root:        ${OUT_ROOT}"
echo "  out csv:         ${OUT_CSV}"
echo "  out slice csv:   ${OUT_SLICE_CSV}"
echo "  all plots pdf:   ${ALL_PLOTS_PDF}"
echo "============================================================================"

BUILD_DIR="$(mktemp -d "${OUT_DIR}/.build_xsec.XXXXXX")"
trap 'rm -rf "${BUILD_DIR}"' EXIT
XSEC_BIN="${BUILD_DIR}/excl_xsec_no_simc_model"

echo "[build] Compiling xsec executable"
"${CXX_CMD}" "${XSEC_SRC}" $(root-config --cflags --libs) -O2 -std=c++17 -o "${XSEC_BIN}"

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

echo "[run] Running xsec extraction"
"${xsec_cmd[@]}"

echo "============================================================================"
echo "Xsec pipeline complete"
echo "============================================================================"
