#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

ANALYSIS_DRIVER="${REPO_ROOT}/src/analysis/run_parallel_nps_analysis_main.sh"
SIM_GENERATOR="${REPO_ROOT}/scripts/generate_simc_infiles.py"
SIM_CHAIN_WRAPPER="${REPO_ROOT}/scripts/run_simulation_chain.py"
DIAG_INDEX_GEN="${REPO_ROOT}/scripts/generate_diagnostics_index.py"
DIAG_COVERAGE_GEN="${REPO_ROOT}/scripts/generate_diagnostics_coverage.py"
DIAG_REPORT_GEN="${REPO_ROOT}/scripts/generate_modular_diagnostics_reports.py"
MASTER_CONFIG_DEFAULT="${REPO_ROOT}/config/nps_dvcs_all_kins_main.csv"
SIM_CONFIG_DEFAULT="${REPO_ROOT}/config/nps_simulation_kinematics.csv"
SIMC_TEMPLATE_DEFAULT="/u/group/nps/singhav/simc_gfortran_updated/infiles/nps_excl_pi0_x60_4b.inp"
SIMC_OUTDIR_DEFAULT="${REPO_ROOT}/config/simc_infiles"
OUTPUT_BASE_DEFAULT="${REPO_ROOT}/output"

CONFIG_CSV="${MASTER_CONFIG_DEFAULT}"
SIM_CONFIG_CSV="${SIM_CONFIG_DEFAULT}"
OUTPUT_BASE="${OUTPUT_BASE_DEFAULT}"
GENERATE_SIM_CONFIG="yes"
GENERATE_SIMC_INFILES="no"
SIMC_TEMPLATE="${SIMC_TEMPLATE_DEFAULT}"
SIMC_OUTDIR="${SIMC_OUTDIR_DEFAULT}"
SIMC_INFILE_PREFIX="nps_excl_pi0_"
SIMC_KINS=()
RUN_SIM_CHAIN="no"
SIM_CHAIN_RUN_SIMC="no"
SIM_CHAIN_RUN_GEANT4="no"
SIM_CHAIN_SIMC_CMD_TEMPLATE=""
SIM_CHAIN_GEANT4_CMD_TEMPLATE=""
SIM_CHAIN_MANIFEST=""
SIM_CHAIN_KINS=()
BUILD_DIAGNOSTICS_INDEX="no"
DIAG_INDEX_OUT=""
BUILD_DIAGNOSTICS_COVERAGE="no"
DIAG_COVERAGE_OUT=""
BUILD_HMS_DIAGNOSTICS="no"
BUILD_NPS_DIAGNOSTICS="no"
BUILD_EXPERIMENT_DIAGNOSTICS="no"
BUILD_DIAGNOSTICS_REPORTS="no"
DIAG_HMS_OUT=""
DIAG_NPS_OUT=""
DIAG_EXPERIMENT_OUT=""
SELECTION_RULE=""

ANALYSIS_ARGS=()

usage() {
  cat <<EOF
Usage: $(basename "$0") [options] [-- <extra args for analysis driver>]

High-level options:
  --config <path>            Master config CSV (default: ${MASTER_CONFIG_DEFAULT})
  --sim-config <path>        Simulation kinematics CSV output path (default: ${SIM_CONFIG_DEFAULT})
  --output-base <path>       Output base path (default: ${OUTPUT_BASE_DEFAULT})
  --generate-sim-config      Generate/update simulation config CSV before running analysis (default)
  --no-generate-sim-config   Skip simulation config generation
  --generate-simc-infiles    Generate infiles from reviewed CSV; skips CSV regeneration
  --no-generate-simc-infiles Skip SIMC infile generation (default)
  --simc-template <path>     Template SIMC infile used for generation
  --simc-outdir <path>       Output directory for generated SIMC infiles
  --simc-file-prefix <text>  Output filename prefix (default: nps_excl_pi0_)
  --simc-kin <Kin_old>       Restrict SIMC infile generation to selected Kin_old (repeatable)
  --run-sim-chain            Run SIMC/Geant4 wrapper stage from canonical simulation CSV
  --sim-chain-run-simc       Enable SIMC command template execution in wrapper
  --sim-chain-run-geant4     Enable Geant4 command template execution in wrapper
  --sim-chain-simc-cmd <cmd>
                              SIMC command template (must include {simc_infile}, {simc_output_file}, {kin_old})
  --sim-chain-geant4-cmd <cmd>
                              Geant4 command template (must include {kin_old}, {geant4_output_dir})
  --sim-chain-kin <Kin_old>  Restrict wrapper stage to selected Kin_old (repeatable)
  --sim-chain-manifest <path> Override wrapper manifest CSV output path
  --build-diagnostics-index  Generate diagnostics artifact index CSV after analysis
  --diagnostics-index-out <path>
                              Override diagnostics index CSV output path
  --build-diagnostics-coverage
                              Generate HMS/NPS diagnostics coverage CSV after analysis
  --diagnostics-coverage-out <path>
                              Override diagnostics coverage CSV output path
  --build-diagnostics-reports
                              Generate modular diagnostics reports (HMS, NPS, experiment)
  --build-hms-diagnostics     Generate HMS diagnostics report after analysis
  --build-nps-diagnostics     Generate NPS diagnostics report after analysis
  --build-experiment-diagnostics
                              Generate experiment-wide diagnostics report after analysis
  --diagnostics-hms-out <path>
                              Override HMS diagnostics report CSV output path
  --diagnostics-nps-out <path>
                              Override NPS diagnostics report CSV output path
  --diagnostics-experiment-out <path>
                              Override experiment diagnostics report CSV output path
  --selection-rule <text>    Provenance text stored in simulation config CSV
  --help                     Show help

Any argument not consumed above is forwarded to:
  ${ANALYSIS_DRIVER}

Examples:
  $(basename "$0") --kin KinC_x60_4b --source waveform --jobs 4 --gevnum-cut no
  $(basename "$0") --no-generate-sim-config -- --all-kins --run-smearing --run-xsec
EOF
}

if [[ ! -f "${ANALYSIS_DRIVER}" ]]; then
  echo "Analysis driver not found: ${ANALYSIS_DRIVER}" >&2
  exit 1
fi
if [[ ! -f "${SIM_GENERATOR}" ]]; then
  echo "Simulation input generator not found: ${SIM_GENERATOR}" >&2
  exit 1
fi
if [[ ! -f "${SIM_CHAIN_WRAPPER}" ]]; then
  echo "Simulation chain wrapper not found: ${SIM_CHAIN_WRAPPER}" >&2
  exit 1
fi
if [[ ! -f "${DIAG_INDEX_GEN}" ]]; then
  echo "Diagnostics index generator not found: ${DIAG_INDEX_GEN}" >&2
  exit 1
fi
if [[ ! -f "${DIAG_COVERAGE_GEN}" ]]; then
  echo "Diagnostics coverage generator not found: ${DIAG_COVERAGE_GEN}" >&2
  exit 1
fi
if [[ ! -f "${DIAG_REPORT_GEN}" ]]; then
  echo "Modular diagnostics report generator not found: ${DIAG_REPORT_GEN}" >&2
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_CSV="$2"
      shift 2
      ;;
    --sim-config)
      SIM_CONFIG_CSV="$2"
      shift 2
      ;;
    --output-base)
      OUTPUT_BASE="$2"
      shift 2
      ;;
    --generate-sim-config)
      GENERATE_SIM_CONFIG="yes"
      shift
      ;;
    --no-generate-sim-config)
      GENERATE_SIM_CONFIG="no"
      shift
      ;;
    --generate-simc-infiles)
      GENERATE_SIMC_INFILES="yes"
      GENERATE_SIM_CONFIG="no"
      shift
      ;;
    --no-generate-simc-infiles)
      GENERATE_SIMC_INFILES="no"
      shift
      ;;
    --simc-template)
      SIMC_TEMPLATE="$2"
      shift 2
      ;;
    --simc-outdir)
      SIMC_OUTDIR="$2"
      shift 2
      ;;
    --simc-file-prefix)
      SIMC_INFILE_PREFIX="$2"
      shift 2
      ;;
    --simc-kin)
      SIMC_KINS+=("$2")
      shift 2
      ;;
    --run-sim-chain)
      RUN_SIM_CHAIN="yes"
      shift
      ;;
    --sim-chain-run-simc)
      SIM_CHAIN_RUN_SIMC="yes"
      shift
      ;;
    --sim-chain-run-geant4)
      SIM_CHAIN_RUN_GEANT4="yes"
      shift
      ;;
    --sim-chain-simc-cmd)
      SIM_CHAIN_SIMC_CMD_TEMPLATE="$2"
      shift 2
      ;;
    --sim-chain-geant4-cmd)
      SIM_CHAIN_GEANT4_CMD_TEMPLATE="$2"
      shift 2
      ;;
    --sim-chain-kin)
      SIM_CHAIN_KINS+=("$2")
      shift 2
      ;;
    --sim-chain-manifest)
      SIM_CHAIN_MANIFEST="$2"
      shift 2
      ;;
    --build-diagnostics-index)
      BUILD_DIAGNOSTICS_INDEX="yes"
      shift
      ;;
    --diagnostics-index-out)
      DIAG_INDEX_OUT="$2"
      shift 2
      ;;
    --build-diagnostics-coverage)
      BUILD_DIAGNOSTICS_COVERAGE="yes"
      shift
      ;;
    --diagnostics-coverage-out)
      DIAG_COVERAGE_OUT="$2"
      shift 2
      ;;
    --build-diagnostics-reports)
      BUILD_DIAGNOSTICS_REPORTS="yes"
      shift
      ;;
    --build-hms-diagnostics)
      BUILD_HMS_DIAGNOSTICS="yes"
      shift
      ;;
    --build-nps-diagnostics)
      BUILD_NPS_DIAGNOSTICS="yes"
      shift
      ;;
    --build-experiment-diagnostics)
      BUILD_EXPERIMENT_DIAGNOSTICS="yes"
      shift
      ;;
    --diagnostics-hms-out)
      DIAG_HMS_OUT="$2"
      shift 2
      ;;
    --diagnostics-nps-out)
      DIAG_NPS_OUT="$2"
      shift 2
      ;;
    --diagnostics-experiment-out)
      DIAG_EXPERIMENT_OUT="$2"
      shift 2
      ;;
    --selection-rule)
      SELECTION_RULE="$2"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    --)
      shift
      while [[ $# -gt 0 ]]; do
        ANALYSIS_ARGS+=("$1")
        shift
      done
      ;;
    *)
      ANALYSIS_ARGS+=("$1")
      shift
      ;;
  esac
done

if [[ ! -f "${CONFIG_CSV}" ]]; then
  echo "Config CSV not found: ${CONFIG_CSV}" >&2
  exit 1
fi

if [[ "${GENERATE_SIMC_INFILES}" == "yes" ]]; then
  if [[ ! -f "${SIMC_TEMPLATE}" ]]; then
    echo "SIMC template not found: ${SIMC_TEMPLATE}" >&2
    exit 1
  fi

  if [[ ${#SIMC_KINS[@]} -eq 0 ]]; then
    # If not explicitly set, mirror forwarded analysis kinematic selection.
    i=0
    while [[ $i -lt ${#ANALYSIS_ARGS[@]} ]]; do
      arg="${ANALYSIS_ARGS[$i]}"
      if [[ "${arg}" == "--kin" ]] && [[ $((i + 1)) -lt ${#ANALYSIS_ARGS[@]} ]]; then
        SIMC_KINS+=("${ANALYSIS_ARGS[$((i + 1))]}")
        i=$((i + 2))
        continue
      fi
      if [[ "${arg}" == --kin=* ]]; then
        SIMC_KINS+=("${arg#--kin=}")
      fi
      i=$((i + 1))
    done
  fi

  SIMC_KIN_ARGS=()
  for kin in "${SIMC_KINS[@]}"; do
    [[ -n "${kin}" ]] || continue
    SIMC_KIN_ARGS+=(--kin "${kin}")
  done

fi

if [[ "${GENERATE_SIM_CONFIG}" == "yes" ]]; then
  echo "[pipeline] Generating simulation kinematics CSV"
  python3 "${SIM_GENERATOR}" \
    --config "${CONFIG_CSV}" \
    --sim-config "${SIM_CONFIG_CSV}" \
    --source-simc-infile "${SIMC_TEMPLATE}" \
    --selection-rule "${SELECTION_RULE}"
fi

if [[ "${GENERATE_SIMC_INFILES}" == "yes" ]]; then
  echo "[pipeline] Generating SIMC infiles from reviewed CSV"
  python3 "${SIM_GENERATOR}" \
    --gen_infile \
    --sim-config "${SIM_CONFIG_CSV}" \
    --template "${SIMC_TEMPLATE}" \
    --out-dir "${SIMC_OUTDIR}" \
    --channel all \
    --file-prefix "${SIMC_INFILE_PREFIX}" \
    "${SIMC_KIN_ARGS[@]}"
fi

if [[ "${RUN_SIM_CHAIN}" == "yes" ]]; then
  if [[ "${SIM_CHAIN_RUN_SIMC}" != "yes" && "${SIM_CHAIN_RUN_GEANT4}" != "yes" ]]; then
    echo "--run-sim-chain requires --sim-chain-run-simc and/or --sim-chain-run-geant4" >&2
    exit 1
  fi

  if [[ ${#SIM_CHAIN_KINS[@]} -eq 0 ]]; then
    # Mirror forwarded analysis --kin selection for wrapper stage by default.
    i=0
    while [[ $i -lt ${#ANALYSIS_ARGS[@]} ]]; do
      arg="${ANALYSIS_ARGS[$i]}"
      if [[ "${arg}" == "--kin" ]] && [[ $((i + 1)) -lt ${#ANALYSIS_ARGS[@]} ]]; then
        SIM_CHAIN_KINS+=("${ANALYSIS_ARGS[$((i + 1))]}")
        i=$((i + 2))
        continue
      fi
      if [[ "${arg}" == --kin=* ]]; then
        SIM_CHAIN_KINS+=("${arg#--kin=}")
      fi
      i=$((i + 1))
    done
  fi

  SIM_CHAIN_KIN_ARGS=()
  for kin in "${SIM_CHAIN_KINS[@]}"; do
    [[ -n "${kin}" ]] || continue
    SIM_CHAIN_KIN_ARGS+=(--kin "${kin}")
  done

  SIM_CHAIN_STAGE_ARGS=(
    --sim-config "${SIM_CONFIG_CSV}"
    --simc-infile-dir "${SIMC_OUTDIR}"
    --simc-file-prefix "${SIMC_INFILE_PREFIX}"
    --output-base "${OUTPUT_BASE}"
    "${SIM_CHAIN_KIN_ARGS[@]}"
  )

  if [[ "${SIM_CHAIN_RUN_SIMC}" == "yes" ]]; then
    SIM_CHAIN_STAGE_ARGS+=(--run-simc --simc-cmd-template "${SIM_CHAIN_SIMC_CMD_TEMPLATE}")
  fi
  if [[ "${SIM_CHAIN_RUN_GEANT4}" == "yes" ]]; then
    SIM_CHAIN_STAGE_ARGS+=(--run-geant4 --geant4-cmd-template "${SIM_CHAIN_GEANT4_CMD_TEMPLATE}")
  fi
  if [[ -n "${SIM_CHAIN_MANIFEST}" ]]; then
    SIM_CHAIN_STAGE_ARGS+=(--manifest "${SIM_CHAIN_MANIFEST}")
  fi

  echo "[pipeline] Running simulation chain wrapper stage"
  python3 "${SIM_CHAIN_WRAPPER}" "${SIM_CHAIN_STAGE_ARGS[@]}"
fi

echo "[pipeline] Running analysis driver"
bash "${ANALYSIS_DRIVER}" \
  --config "${CONFIG_CSV}" \
  --output-base "${OUTPUT_BASE}" \
  "${ANALYSIS_ARGS[@]}"
rc=$?
if [[ ${rc} -ne 0 ]]; then
  echo "[pipeline] Analysis driver failed with exit code ${rc}" >&2
  exit ${rc}
fi

if [[ "${BUILD_DIAGNOSTICS_INDEX}" == "yes" ]]; then
  DIAG_INDEX_ARGS=(--output-base "${OUTPUT_BASE}")
  if [[ -n "${DIAG_INDEX_OUT}" ]]; then
    DIAG_INDEX_ARGS+=(--out-file "${DIAG_INDEX_OUT}")
  fi
  echo "[pipeline] Generating diagnostics artifact index"
  python3 "${DIAG_INDEX_GEN}" "${DIAG_INDEX_ARGS[@]}"
fi

if [[ "${BUILD_DIAGNOSTICS_COVERAGE}" == "yes" ]]; then
  DIAG_COVERAGE_ARGS=(--output-base "${OUTPUT_BASE}")
  if [[ -n "${DIAG_COVERAGE_OUT}" ]]; then
    DIAG_COVERAGE_ARGS+=(--out-file "${DIAG_COVERAGE_OUT}")
  fi
  echo "[pipeline] Generating HMS/NPS diagnostics coverage"
  python3 "${DIAG_COVERAGE_GEN}" "${DIAG_COVERAGE_ARGS[@]}"
fi

# Backward-compatible mapping: old coverage flag implies experiment diagnostics report.
if [[ "${BUILD_DIAGNOSTICS_COVERAGE}" == "yes" ]]; then
  BUILD_EXPERIMENT_DIAGNOSTICS="yes"
  if [[ -n "${DIAG_COVERAGE_OUT}" && -z "${DIAG_EXPERIMENT_OUT}" ]]; then
    DIAG_EXPERIMENT_OUT="${DIAG_COVERAGE_OUT}"
  fi
fi

if [[ "${BUILD_DIAGNOSTICS_REPORTS}" == "yes" ]]; then
  BUILD_HMS_DIAGNOSTICS="yes"
  BUILD_NPS_DIAGNOSTICS="yes"
  BUILD_EXPERIMENT_DIAGNOSTICS="yes"
fi

if [[ "${BUILD_HMS_DIAGNOSTICS}" == "yes" || "${BUILD_NPS_DIAGNOSTICS}" == "yes" || "${BUILD_EXPERIMENT_DIAGNOSTICS}" == "yes" ]]; then
  if [[ "${BUILD_HMS_DIAGNOSTICS}" == "yes" ]]; then
    DIAG_HMS_ARGS=(--report hms)
    if [[ -n "${DIAG_HMS_OUT}" ]]; then
      DIAG_HMS_ARGS+=(--hms-out "${DIAG_HMS_OUT}")
    fi
    echo "[pipeline] Generating HMS diagnostics report"
    python3 "${DIAG_REPORT_GEN}" --output-base "${OUTPUT_BASE}" "${DIAG_HMS_ARGS[@]}"
  fi

  if [[ "${BUILD_NPS_DIAGNOSTICS}" == "yes" ]]; then
    DIAG_NPS_ARGS=(--report nps)
    if [[ -n "${DIAG_NPS_OUT}" ]]; then
      DIAG_NPS_ARGS+=(--nps-out "${DIAG_NPS_OUT}")
    fi
    echo "[pipeline] Generating NPS diagnostics report"
    python3 "${DIAG_REPORT_GEN}" --output-base "${OUTPUT_BASE}" "${DIAG_NPS_ARGS[@]}"
  fi

  if [[ "${BUILD_EXPERIMENT_DIAGNOSTICS}" == "yes" ]]; then
    DIAG_EXPERIMENT_ARGS=(--report experiment)
    if [[ -n "${DIAG_EXPERIMENT_OUT}" ]]; then
      DIAG_EXPERIMENT_ARGS+=(--experiment-out "${DIAG_EXPERIMENT_OUT}")
    fi
    echo "[pipeline] Generating experiment diagnostics report"
    python3 "${DIAG_REPORT_GEN}" --output-base "${OUTPUT_BASE}" "${DIAG_EXPERIMENT_ARGS[@]}"
  fi
fi
