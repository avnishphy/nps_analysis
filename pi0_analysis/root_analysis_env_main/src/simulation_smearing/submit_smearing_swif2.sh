#!/usr/bin/env bash
# Build and optionally submit one NPS smearing job to ifarm through SWIF2.
#
# The default job is preview-only: run_smearing_pipeline.sh writes the input
# comparison PDF, receives EOF at its approval prompt, and stops before fitting.
# Set INPUTS_APPROVED=1 (or pass --inputs-approved) only after reviewing that PDF.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PIPELINE="${SCRIPT_DIR}/run_smearing_pipeline.sh"

ACCOUNT="${ACCOUNT:-hallc}"
PARTITION="${PARTITION:-production}"
WORKFLOW="${WORKFLOW:-nps_pi0_smearing_$(date +%Y%m%d_%H%M%S)}"
JOB_NAME="${JOB_NAME:-${WORKFLOW}_job0}"

JSON_DIR="${JSON_DIR:-/group/nps/${USER}/hcswif/jsons}"
STDOUT_DIR="${STDOUT_DIR:-/farm_out/${USER}/nps_replay_stdout}"
STDERR_DIR="${STDERR_DIR:-/farm_out/${USER}/nps_replay_stderr}"
WRAPPER_DIR="${WRAPPER_DIR:-${SCRIPT_DIR}/farm_jobs}"

CPU_CORES="${CPU_CORES:-1}"
OMP_THREADS="${OMP_THREADS:-${CPU_CORES}}"
RAM_BYTES="${RAM_BYTES:-64000000000}"
TIME_SECS="${TIME_SECS:-172800}"

KIN="${KIN:-${NPS_SMEAR_KIN:-x60_4b}}"
TARGET="${TARGET:-${NPS_SMEAR_TARGET:-LH2}}"
ROOT_DIR="${ROOT_DIR:-}"
COMBINED_FILE="${COMBINED_FILE:-}"
SIM_FILE="${SIM_FILE:-}"
SIM_SMEARED_FILE="${SIM_SMEARED_FILE:-}"
OUT_FILE="${OUT_FILE:-}"
SECTION_MAP_FILE="${SECTION_MAP_FILE:-}"
INTERP_FILE="${INTERP_FILE:-}"
INPUT_PREVIEW_PDF="${INPUT_PREVIEW_PDF:-}"

NX="${NX:-11}"
NY="${NY:-16}"
X_MIN="${X_MIN:--24}"
X_MAX="${X_MAX:-28}"
Y_MIN="${Y_MIN:--34}"
Y_MAX="${Y_MAX:-34}"
OVERLAP="${OVERLAP:-0.0}"
NSMEAR="${NSMEAR:-80}"
SMEAR_SEED="${SMEAR_SEED:-42}"

BEAM_ENERGY_GEV="${BEAM_ENERGY_GEV:-${NPS_EBEAM:-}}"
Z_NPS_CM="${Z_NPS_CM:-${NPS_Z_NPS_CM:-}}"
THETA_NPS_DEG="${THETA_NPS_DEG:-${NPS_THETA_NPS_DEG:-}}"
NPS_ANGLE_DEG="${NPS_ANGLE_DEG:-${NPS_NPS_THETA_DEG:-}}"
MODE="${MODE:-auto}"
ACCEPTANCE_CONFIG="${ACCEPTANCE_CONFIG:-}"
SIM_EXCL_INPUT="${SIM_EXCL_INPUT:-}"
SIM_SIDIS_INPUT="${SIM_SIDIS_INPUT:-}"
SIM_DELTA_INPUT="${SIM_DELTA_INPUT:-}"

# Safety controls. JSON_ONLY avoids SWIF2 state changes; INPUTS_APPROVED controls
# whether the non-interactive job may answer yes to the fitter's approval prompt.
JSON_ONLY="${JSON_ONLY:-0}"
INPUTS_APPROVED="${INPUTS_APPROVED:-0}"

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --combined-file <path>       Combined data ROOT file
  --beam-energy-gev <value>    Beam energy in GeV
  --z-nps-cm <value>           Target-to-NPS distance in cm
  --theta-nps-deg <value>      Signed NPS-to-Hall angle, or use --nps-angle-deg

Physics and output options:
  --kin <name>                 Kinematic name (default: x60_4b)
  --target <name>              Target token (default: LH2)
  --root-dir <path>            Pipeline ROOT/output directory
  --sim-file <path>            Unsmeared SIMC ROOT file
  --sim-smeared-file <path>    Smeared SIMC ROOT output
  --out-file <path>            Smearing-fit ROOT output
  --section-map-file <path>    Section-map CSV output
  --interp-file <path>         Interpolated-map ROOT output
  --input-preview-pdf <path>   Input-comparison PDF
  --nps-angle-deg <value>      Physical NPS angle (alternative to --theta-nps-deg)
  --nx/--ny <int>              Section counts (defaults: 11 and 16)
  --x-min/--x-max <value>      X bounds (defaults: -24 and 28)
  --y-min/--y-max <value>      Y bounds (defaults: -34 and 34)
  --overlap <value>            Section overlap fraction (default: 0.0)
  --nsmear <int>               Smearing iterations/event (default: 80)
  --smear-seed <int>           SIMC smearing RNG seed (default: 42)
  --mode <value>               SIMC mode hint (default: auto)
  --acceptance-config <path>   Acceptance-cuts configuration
  --sim-excl-input <path>      Exclusive SIMC input
  --sim-sidis-input <path>     SIDIS SIMC input
  --sim-delta-input <path>     Delta-pi0 SIMC input

Submission controls:
  --workflow <name>            SWIF2 workflow name
  --job-name <name>            SWIF2 job name
  --inputs-approved            Answer yes to the fitter prompt and run the fit
  --preview-only               Generate preview PDF and stop before fitting (default)
  --json-only                  Generate wrapper/JSON without importing or running
  --help                       Show this message

Resources and farm paths can be overridden with ACCOUNT, PARTITION, CPU_CORES,
OMP_THREADS, RAM_BYTES, TIME_SECS, JSON_DIR, STDOUT_DIR, STDERR_DIR, and
WRAPPER_DIR. All physics options also accept their uppercase environment names.

Examples:
  # First submission: generate the PDF only.
  $(basename "$0") --kin KinC_x36_5 --combined-file /path/data.root \\
    --beam-energy-gev 10.538 --z-nps-cm 407 --nps-angle-deg 12.2

  # After reviewing the PDF, submit the optimization explicitly.
  $(basename "$0") --kin KinC_x36_5 --combined-file /path/data.root \\
    --beam-energy-gev 10.538 --z-nps-cm 407 --nps-angle-deg 12.2 \\
    --inputs-approved
EOF
}

need_value() {
  if [[ $# -lt 2 || -z "${2}" ]]; then
    echo "[ERROR] $1 requires a value." >&2
    exit 2
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --kin) need_value "$@"; KIN="$2"; shift 2 ;;
    --target) need_value "$@"; TARGET="$2"; shift 2 ;;
    --root-dir) need_value "$@"; ROOT_DIR="$2"; shift 2 ;;
    --combined-file) need_value "$@"; COMBINED_FILE="$2"; shift 2 ;;
    --sim-file) need_value "$@"; SIM_FILE="$2"; shift 2 ;;
    --sim-smeared-file) need_value "$@"; SIM_SMEARED_FILE="$2"; shift 2 ;;
    --out-file) need_value "$@"; OUT_FILE="$2"; shift 2 ;;
    --section-map-file) need_value "$@"; SECTION_MAP_FILE="$2"; shift 2 ;;
    --interp-file) need_value "$@"; INTERP_FILE="$2"; shift 2 ;;
    --input-preview-pdf) need_value "$@"; INPUT_PREVIEW_PDF="$2"; shift 2 ;;
    --nx) need_value "$@"; NX="$2"; shift 2 ;;
    --ny) need_value "$@"; NY="$2"; shift 2 ;;
    --x-min) need_value "$@"; X_MIN="$2"; shift 2 ;;
    --x-max) need_value "$@"; X_MAX="$2"; shift 2 ;;
    --y-min) need_value "$@"; Y_MIN="$2"; shift 2 ;;
    --y-max) need_value "$@"; Y_MAX="$2"; shift 2 ;;
    --overlap) need_value "$@"; OVERLAP="$2"; shift 2 ;;
    --nsmear) need_value "$@"; NSMEAR="$2"; shift 2 ;;
    --smear-seed) need_value "$@"; SMEAR_SEED="$2"; shift 2 ;;
    --beam-energy-gev) need_value "$@"; BEAM_ENERGY_GEV="$2"; shift 2 ;;
    --z-nps-cm) need_value "$@"; Z_NPS_CM="$2"; shift 2 ;;
    --theta-nps-deg) need_value "$@"; THETA_NPS_DEG="$2"; shift 2 ;;
    --nps-angle-deg) need_value "$@"; NPS_ANGLE_DEG="$2"; shift 2 ;;
    --mode) need_value "$@"; MODE="$2"; shift 2 ;;
    --acceptance-config) need_value "$@"; ACCEPTANCE_CONFIG="$2"; shift 2 ;;
    --sim-excl-input) need_value "$@"; SIM_EXCL_INPUT="$2"; shift 2 ;;
    --sim-sidis-input) need_value "$@"; SIM_SIDIS_INPUT="$2"; shift 2 ;;
    --sim-delta-input) need_value "$@"; SIM_DELTA_INPUT="$2"; shift 2 ;;
    --workflow) need_value "$@"; WORKFLOW="$2"; shift 2 ;;
    --job-name) need_value "$@"; JOB_NAME="$2"; shift 2 ;;
    --inputs-approved) INPUTS_APPROVED=1; shift ;;
    --preview-only) INPUTS_APPROVED=0; shift ;;
    --json-only) JSON_ONLY=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "[ERROR] Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

require_file() {
  local path="$1"
  local label="$2"
  if [[ ! -f "$path" ]]; then
    echo "[ERROR] Missing ${label}: ${path}" >&2
    exit 1
  fi
}

require_command() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "[ERROR] Command not found: $1" >&2
    exit 1
  fi
}

require_uint() {
  if ! [[ "$1" =~ ^[0-9]+$ ]]; then
    echo "[ERROR] $2 must be a non-negative integer (got: $1)" >&2
    exit 1
  fi
}

if [[ -z "$COMBINED_FILE" || -z "$BEAM_ENERGY_GEV" || -z "$Z_NPS_CM" ]]; then
  echo "[ERROR] --combined-file, --beam-energy-gev, and --z-nps-cm are required." >&2
  exit 2
fi
if [[ -n "$THETA_NPS_DEG" && -n "$NPS_ANGLE_DEG" ]]; then
  echo "[ERROR] Use either --theta-nps-deg or --nps-angle-deg, not both." >&2
  exit 2
fi
if [[ -z "$THETA_NPS_DEG" && -z "$NPS_ANGLE_DEG" ]]; then
  echo "[ERROR] --theta-nps-deg or --nps-angle-deg is required." >&2
  exit 2
fi
if [[ "$INPUTS_APPROVED" != 0 && "$INPUTS_APPROVED" != 1 ]]; then
  echo "[ERROR] INPUTS_APPROVED must be 0 or 1." >&2
  exit 2
fi
if [[ "$JSON_ONLY" != 0 && "$JSON_ONLY" != 1 ]]; then
  echo "[ERROR] JSON_ONLY must be 0 or 1." >&2
  exit 2
fi

for entry in "CPU_CORES:$CPU_CORES" "OMP_THREADS:$OMP_THREADS" \
             "RAM_BYTES:$RAM_BYTES" "TIME_SECS:$TIME_SECS" \
             "NX:$NX" "NY:$NY" "NSMEAR:$NSMEAR" "SMEAR_SEED:$SMEAR_SEED"; do
  require_uint "${entry#*:}" "${entry%%:*}"
done

require_file "$PIPELINE" "smearing pipeline"
require_file "$COMBINED_FILE" "combined data ROOT file"
require_command python3

WRAPPER="${WRAPPER_DIR}/run_smearing_${WORKFLOW}.sh"
JSON_FILE="${JSON_DIR}/${WORKFLOW}.json"

pipeline_args=(
  "$PIPELINE"
  --kin "$KIN"
  --target "$TARGET"
  --combined-file "$COMBINED_FILE"
  --nx "$NX" --ny "$NY"
  --x-min "$X_MIN" --x-max "$X_MAX"
  --y-min "$Y_MIN" --y-max "$Y_MAX"
  --overlap "$OVERLAP"
  --nsmear "$NSMEAR"
  --smear-seed "$SMEAR_SEED"
  --beam-energy-gev "$BEAM_ENERGY_GEV"
  --z-nps-cm "$Z_NPS_CM"
  --mode "$MODE"
)

if [[ -n "$THETA_NPS_DEG" ]]; then
  pipeline_args+=(--theta-nps-deg "$THETA_NPS_DEG")
else
  pipeline_args+=(--nps-angle-deg "$NPS_ANGLE_DEG")
fi

optional_args=(
  ROOT_DIR --root-dir
  SIM_FILE --sim-file
  SIM_SMEARED_FILE --sim-smeared-file
  OUT_FILE --out-file
  SECTION_MAP_FILE --section-map-file
  INTERP_FILE --interp-file
  INPUT_PREVIEW_PDF --input-preview-pdf
  ACCEPTANCE_CONFIG --acceptance-config
  SIM_EXCL_INPUT --sim-excl-input
  SIM_SIDIS_INPUT --sim-sidis-input
  SIM_DELTA_INPUT --sim-delta-input
)
for ((i = 0; i < ${#optional_args[@]}; i += 2)); do
  var_name="${optional_args[i]}"
  option_name="${optional_args[i + 1]}"
  if [[ -n "${!var_name}" ]]; then
    pipeline_args+=("$option_name" "${!var_name}")
  fi
done

printf -v pipeline_command '%q ' "${pipeline_args[@]}"

mkdir -p "$JSON_DIR" "$STDOUT_DIR" "$STDERR_DIR" "$WRAPPER_DIR"

cat > "$WRAPPER" <<EOF
#!/usr/bin/env bash
set -euo pipefail

# Reload this wrapper once inside the Hall C software environment.
if [[ "\${NPS_FARM_ENV_LOADED:-0}" != 1 ]]; then
  export NPS_FARM_ENV_LOADED=1
  exec csh -c 'if ( -f /usr/share/Modules/init/csh ) source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; bash "$WRAPPER" 0'
fi

export OMP_NUM_THREADS="$OMP_THREADS"
cd "$REPO_ROOT"

if [[ "$INPUTS_APPROVED" == 1 ]]; then
  printf 'y\\n' | $pipeline_command
else
  $pipeline_command </dev/null
fi
EOF
chmod +x "$WRAPPER"

python3 - "$JSON_FILE" <<PY
import json
import sys

json_file = sys.argv[1]
workflow = ${WORKFLOW@Q}
job_name = ${JOB_NAME@Q}
wrapper = ${WRAPPER@Q}
stdout = ${STDOUT_DIR@Q} + "/" + job_name + ".out"
stderr = ${STDERR_DIR@Q} + "/" + job_name + ".err"

payload = {
    "name": workflow,
    "jobs": [
        {
            "name": job_name,
            "account": ${ACCOUNT@Q},
            "partition": ${PARTITION@Q},
            "command": [wrapper + " 0"],
            "cpu_cores": int(${CPU_CORES@Q}),
            "ram_bytes": int(${RAM_BYTES@Q}),
            "time_secs": int(${TIME_SECS@Q}),
            "stdout": stdout,
            "stderr": stderr,
        }
    ],
}

with open(json_file, "w", encoding="utf-8") as stream:
    json.dump(payload, stream, indent=2)
    stream.write("\\n")
PY

echo "Workflow:        $WORKFLOW"
echo "Job:             $JOB_NAME"
echo "Mode:            $([[ "$INPUTS_APPROVED" == 1 ]] && echo optimization-approved || echo preview-only)"
echo "Wrapper:         $WRAPPER"
echo "JSON:            $JSON_FILE"
echo "Stdout:          $STDOUT_DIR/$JOB_NAME.out"
echo "Stderr:          $STDERR_DIR/$JOB_NAME.err"
echo "CPU/RAM/time:    $CPU_CORES / $RAM_BYTES / $TIME_SECS"

if [[ "$JSON_ONLY" == 1 ]]; then
  echo "JSON_ONLY=1: wrapper and JSON generated; SWIF2 was not changed."
  exit 0
fi

require_command swif2
swif2 import -file "$JSON_FILE"
swif2 run "$WORKFLOW"

echo "Submitted SWIF2 workflow: $WORKFLOW"
echo "Check status with: swif2 status $WORKFLOW"
