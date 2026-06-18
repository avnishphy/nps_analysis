#!/bin/bash
# Submit the pi0 simulation smearing pipeline to ifarm with SWIF2.
#
# Typical use:
#   ./submit_job_smearing.sh
#
# Useful overrides:
#   NX=11 NY=16 NSMEAR=80 ./submit_job_smearing.sh
#   OUT_DIR=/volatile/hallc/nps/$USER/pi0_smearing/test JSON_ONLY=1 ./submit_job_smearing.sh

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

ACCOUNT="${ACCOUNT:-hallc}"
PARTITION="${PARTITION:-production}"
WORKFLOW="${WORKFLOW:-nps_pi0_smearing_$(date +%Y%m%d_%H%M%S)}"
JOB_NAME="${JOB_NAME:-${WORKFLOW}_job0}"

JSON_DIR="${JSON_DIR:-/group/nps/$USER/hcswif/jsons}"
STDOUT_DIR="${STDOUT_DIR:-/farm_out/$USER/nps_replay_stdout}"
STDERR_DIR="${STDERR_DIR:-/farm_out/$USER/nps_replay_stderr}"
WRAPPER_DIR="${WRAPPER_DIR:-$PROJECT_ROOT/farm_jobs}"
WRAPPER="$WRAPPER_DIR/run_smearing_${WORKFLOW}.sh"
JSON_FILE="$JSON_DIR/$WORKFLOW.json"

# Keep the farm request serial by default. Override both if you intentionally
# want to run the fitter with more OpenMP threads.
CPU_CORES="${CPU_CORES:-1}"
OMP_THREADS="${OMP_THREADS:-$CPU_CORES}"

# Resource defaults can be overridden without changing this script.
RAM_BYTES="${RAM_BYTES:-64000000000}"
TIME_SECS="${TIME_SECS:-172800}"

DATA_FILE="${DATA_FILE:-$PROJECT_ROOT/output/plots/x60_4b/production_wfpi0/combined_branches_LH2_wfpi0.root}"
SIM_FILE="${SIM_FILE:-$PROJECT_ROOT/output/plots/x60_4b/production_wfpi0/simc_pi0_analysis_output.root}"
OUT_DIR="${OUT_DIR:-/volatile/hallc/nps/$USER/pi0_smearing/$WORKFLOW}"
RUN_TAG="${RUN_TAG:-farm_${WORKFLOW}}"

NX="${NX:-7}"
NY="${NY:-7}"
X_MIN="${X_MIN:--24}"
X_MAX="${X_MAX:-28}"
Y_MIN="${Y_MIN:--34}"
Y_MAX="${Y_MAX:-34}"
OVERLAP_FRAC="${OVERLAP_FRAC:-0.0}"
NSMEAR="${NSMEAR:-80}"

RUN_COMPARISON="${RUN_COMPARISON:-1}"
COMPARISON_MAX_ENTRIES="${COMPARISON_MAX_ENTRIES:-}"
NPS_SMEARING_SWEEP_ACCEPTANCE="${NPS_SMEARING_SWEEP_ACCEPTANCE:-jacobi_global_accept_rollback}"

# Set JSON_ONLY=1 to generate the wrapper/JSON without importing or running.
JSON_ONLY="${JSON_ONLY:-0}"

require_file() {
  local path="$1"
  local label="$2"
  if [[ ! -f "$path" ]]; then
    echo "[ERROR] Missing $label: $path" >&2
    exit 1
  fi
}

require_command() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "[ERROR] $cmd not found in PATH" >&2
    exit 1
  fi
}

require_int() {
  local value="$1"
  local label="$2"
  if ! [[ "$value" =~ ^[0-9]+$ ]]; then
    echo "[ERROR] $label must be a non-negative integer: $value" >&2
    exit 1
  fi
}

require_file "$DATA_FILE" "data ROOT file"
require_file "$SIM_FILE" "SIMC ROOT file"
require_file "$PROJECT_ROOT/run_smearing_pipeline.sh" "smearing pipeline"
require_command python3
require_int "$CPU_CORES" "CPU_CORES"
require_int "$OMP_THREADS" "OMP_THREADS"
require_int "$RAM_BYTES" "RAM_BYTES"
require_int "$TIME_SECS" "TIME_SECS"

mkdir -p "$JSON_DIR" "$STDOUT_DIR" "$STDERR_DIR" "$WRAPPER_DIR" "$OUT_DIR"

cat > "$WRAPPER" <<EOF
#!/bin/bash
set -euo pipefail

# Keep a dummy positional argument for compatibility with command-style jobs.
: "\${1:-0}"

export OMP_NUM_THREADS="$OMP_THREADS"
export DATA_FILE="$DATA_FILE"
export SIM_FILE="$SIM_FILE"
export OUT_DIR="$OUT_DIR"
export RUN_TAG="$RUN_TAG"
export NX="$NX"
export NY="$NY"
export X_MIN="$X_MIN"
export X_MAX="$X_MAX"
export Y_MIN="$Y_MIN"
export Y_MAX="$Y_MAX"
export OVERLAP_FRAC="$OVERLAP_FRAC"
export NSMEAR="$NSMEAR"
export RUN_COMPARISON="$RUN_COMPARISON"
export COMPARISON_MAX_ENTRIES="$COMPARISON_MAX_ENTRIES"
export NPS_SMEARING_SWEEP_ACCEPTANCE="$NPS_SMEARING_SWEEP_ACCEPTANCE"

cd "$PROJECT_ROOT"
exec csh -c 'if ( -f /usr/share/Modules/init/csh ) source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; bash ./run_smearing_pipeline.sh'
EOF
chmod +x "$WRAPPER"

python3 - "$JSON_FILE" <<PY
import json
import sys

json_file = sys.argv[1]
workflow = ${WORKFLOW@Q}
job_name = ${JOB_NAME@Q}
wrapper = ${WRAPPER@Q}
account = ${ACCOUNT@Q}
partition = ${PARTITION@Q}
stdout = ${STDOUT_DIR@Q} + "/" + job_name + ".out"
stderr = ${STDERR_DIR@Q} + "/" + job_name + ".err"

payload = {
    "name": workflow,
    "jobs": [
        {
            "name": job_name,
            "account": account,
            "partition": partition,
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

echo "Workflow: $WORKFLOW"
echo "Job:      $JOB_NAME"
echo "Wrapper:  $WRAPPER"
echo "JSON:     $JSON_FILE"
echo "OUT_DIR:  $OUT_DIR"
echo "Stdout:   $STDOUT_DIR/$JOB_NAME.out"
echo "Stderr:   $STDERR_DIR/$JOB_NAME.err"
echo "CPU:      $CPU_CORES"
echo "RAM:      $RAM_BYTES"
echo "TIME:     $TIME_SECS"

if [[ "$JSON_ONLY" == "1" ]]; then
  echo "JSON_ONLY=1, not importing/running SWIF2 workflow."
  exit 0
fi

require_command swif2

swif2 import -file "$JSON_FILE"
swif2 run "$WORKFLOW"

echo "Submitted SWIF2 workflow: $WORKFLOW"
echo "Status: swif2 status $WORKFLOW"
