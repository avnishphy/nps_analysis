#!/usr/bin/env bash
# Run under the sourced Hall C ROOT environment. Creates a fresh campaign.
set -euo pipefail
cd "$(dirname "$0")/.."
PYTHON_CMD="/group/nps/singhav/software/python/bin/python"
export NPS_PYTHON_CMD="${PYTHON_CMD}"
[[ -x "${PYTHON_CMD}" ]] || { echo "Python interpreter is not executable: ${PYTHON_CMD}" >&2; exit 1; }

# Pipeline mode: generate a fresh central/toy campaign from the fit that just
# completed. Every path is explicit so the selected JSON and its derived model
# cache cannot be mixed with a frozen campaign from another bin definition.
if [[ $# -ge 6 && $# -le 8 ]]; then
  OUT="$1"
  FIT_OUTPUT="$2"
  CONFIG="$3"
  DATA_FILE="$4"
  SIM_FILE="$5"
  VERTEX_FILE="$6"
  [[ ! -e "$OUT" ]] || { echo "Refusing existing output: $OUT" >&2; exit 1; }
  command -v root-config >/dev/null
  export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
  TOY_REPLICAS="${7:-500}"
  TOY_JOBS="${8:-0}"
  [[ "${TOY_REPLICAS}" =~ ^[1-9][0-9]*$ ]] || { echo "Toy replica count must be a positive integer" >&2; exit 1; }
  [[ "${TOY_JOBS}" =~ ^[0-9]+$ ]] || { echo "Toy worker count must be zero or a positive integer" >&2; exit 1; }
  mkdir -p "$OUT/before" "$OUT/build"
  cp "$CONFIG" "$OUT/config_snapshot.json"
  cp scripts/preliminary_pi0_xsec.py "$OUT/before/preliminary_pi0_xsec_v7.py"
  "${PYTHON_CMD}" - "$OUT" "$FIT_OUTPUT" "$DATA_FILE" "$SIM_FILE" "$VERTEX_FILE" <<'PY'
from pathlib import Path
import hashlib,json,subprocess,sys,time
out,fit,data,sim,vertex=map(lambda value:Path(value).resolve(),sys.argv[1:])
required=['model_event_cache.csv','model_context.csv','model_parameters.csv','model_reconstructed_yields.csv']
for name in required:
    if not (fit/name).is_file():raise FileNotFoundError(fit/name)
for path in (data,sim,vertex):
    if not path.is_file():raise FileNotFoundError(path)
context=dict(fit_output_dir=str(fit),data_file=str(data),sim_file=str(sim),vertex_file=str(vertex),
             started_epoch=time.time())
(out/'campaign_context.json').write_text(json.dumps(context,indent=2)+'\n')
(out/'before/git_head.txt').write_bytes(subprocess.check_output(['git','rev-parse','HEAD']))
(out/'before/git_status.txt').write_bytes(subprocess.check_output(['git','status','--short']))
(out/'before/input_sha256.json').write_text(json.dumps({str(path):hashlib.sha256(path.read_bytes()).hexdigest()
    for path in [fit/name for name in required]+[out/'config_snapshot.json']},indent=2)+'\n')
PY
  "${PYTHON_CMD}" src/xsec_extract/generate_xsec_config.py "$OUT/config_snapshot.json" "$OUT/build/xsec_config.h"
  g++ -O3 -std=c++17 -shared -fPIC -Isrc/xsec_extract -I"$OUT/build" scripts/pi0_preliminary_bridge.cpp \
    $(root-config --cflags --libs) -lMinuit2 -o "$OUT/build/libprelim.so.tmp"
  mv "$OUT/build/libprelim.so.tmp" "$OUT/build/libprelim.so"
  g++ -O2 -std=c++17 -shared -fPIC src/analysis/nps_stat_bridge.cpp \
    $(root-config --cflags --libs) -lMinuit2 -o "$OUT/build/libnps_stat.so.tmp"
  mv "$OUT/build/libnps_stat.so.tmp" "$OUT/build/libnps_stat.so"
  "${PYTHON_CMD}" scripts/preliminary_pi0_xsec.py central --output "$OUT" > "$OUT/central_final.log" 2>&1
  "${PYTHON_CMD}" scripts/preliminary_pi0_xsec.py toys --output "$OUT" --replicas "$TOY_REPLICAS" \
    --seed 20261007 --jobs "$TOY_JOBS" > "$OUT/toys.log" 2>&1
  "${PYTHON_CMD}" - "$OUT" <<'PY'
from pathlib import Path
import json,sys,time
out=Path(sys.argv[1]);summary=json.loads((out/'toys/summary.json').read_text())
payload=dict(status='complete',accepted_toys=summary['accepted'],finished_epoch=time.time(),
             config_snapshot=str((out/'config_snapshot.json').resolve()),freshly_generated=True)
(out/'generation_summary.json').write_text(json.dumps(payload,indent=2)+'\n')
PY
  exit 0
fi

OUT="${1:-validation/preliminary_model_xsec_reproduction_20261005}"
[[ ! -e "$OUT" ]] || { echo "Refusing existing output: $OUT" >&2; exit 1; }
command -v root-config >/dev/null
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"${PYTHON_CMD}" - "$OUT" <<'PY'
from pathlib import Path
import sys,json,hashlib,subprocess,shutil,time
import platform,numpy,scipy,uproot
out=Path(sys.argv[1]);out.mkdir();(out/'before').mkdir();(out/'build').mkdir()
for cmd,name in [(['git','rev-parse','HEAD'],'git_head.txt'),(['git','status','--short'],'git_status.txt'),(['git','diff','--stat'],'git_diff_stat.txt')]:
 (out/'before'/name).write_bytes(subprocess.check_output(cmd))
files=[p for base in ['src','scripts','config'] for p in Path(base).rglob('*') if p.is_file() and p.suffix in ['.C','.h','.cpp','.py','.sh','.conf','.json','.csv']]
(out/'before/source_sha256.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2))
r=json.loads(Path('validation/model_release_20261004/audit_release/input_audit.json').read_text())['accepted_runs']
(out/'accepted_runs.txt').write_text('\n'.join(map(str,r))+'\n')
shutil.copy2('src/xsec_extract/xsec_config/xsec_config_x36_4.json',out/'config_snapshot.json')
(out/'started_epoch.txt').write_text(str(time.time())+'\n')
(out/'environment.json').write_text(json.dumps(dict(python=sys.version,numpy=numpy.__version__,scipy=scipy.__version__,
 uproot=uproot.__version__,platform=platform.platform(),root=subprocess.check_output(['root-config','--version'],text=True).strip()),indent=2)+'\n')
PY
exec > >(tee "$OUT/driver.log") 2>&1
set -x
bash src/analysis/run_parallel_nps_analysis_main.sh --kin KinC_x36_4 --source waveform \
  --gevnum-cut yes --target LH2 --jobs 3 --types production,Production --no-combine \
  --run $(cat "$OUT/accepted_runs.txt") --output-base "$OUT/production"
"${PYTHON_CMD}" src/analysis/combine_analysis_branches.py --kin KinC_x36_4 --target LH2 \
  --run $(cat "$OUT/accepted_runs.txt") --output-base "$OUT/production" \
  --efficiency-csv output/efficiency_stuff/efficiency_KinC_x36_4.csv --no-analysis-plots
"${PYTHON_CMD}" scripts/validate_pi0_timing_transport.py verify --production-dir "$OUT/production/KinC_x36_4/root" --output "$OUT/production_parity"
"${PYTHON_CMD}" scripts/validate_pi0_timing_transport.py audit --data "$OUT/production/KinC_x36_4/root/combined_branches_LH2.root" --output "$OUT/timing_audit"
"${PYTHON_CMD}" src/xsec_extract/generate_xsec_config.py src/xsec_extract/xsec_config/xsec_config_x36_4.json "$OUT/build/xsec_config.h"
g++ -O3 -std=c++17 -shared -fPIC -Isrc/xsec_extract -I"$OUT/build" scripts/pi0_preliminary_bridge.cpp \
  $(root-config --cflags --libs) -lMinuit2 -o "$OUT/build/libprelim.so.tmp"
mv "$OUT/build/libprelim.so.tmp" "$OUT/build/libprelim.so"
g++ -O2 -std=c++17 -shared -fPIC src/analysis/nps_stat_bridge.cpp \
  $(root-config --cflags --libs) -lMinuit2 -o "$OUT/build/libnps_stat.so.tmp"
mv "$OUT/build/libnps_stat.so.tmp" "$OUT/build/libnps_stat.so"
"${PYTHON_CMD}" scripts/preliminary_pi0_xsec.py central --output "$OUT" > "$OUT/central_final.log" 2>&1
"${PYTHON_CMD}" scripts/validate_preliminary_pi0.py --output "$OUT" > "$OUT/validation.log" 2>&1
"${PYTHON_CMD}" scripts/preliminary_pi0_xsec.py data --output "$OUT" --replicas 1000 --seed 20261005 \
  --stop-after "${PRELIM_DATA_ACCEPTED:-1000}" \
  --stop-reason 'Critical blocker 7: M0 toy pulls are grossly narrower than 0.7; no release.' > "$OUT/data.log" 2>&1 &
DATA_PID=$!
"${PYTHON_CMD}" scripts/preliminary_pi0_xsec.py mc --output "$OUT" --replicas 1000 --seed 20261006 > "$OUT/mc.log" 2>&1 &
MC_PID=$!
"${PYTHON_CMD}" scripts/preliminary_pi0_xsec.py toys --output "$OUT" --replicas 500 --seed 20261007 > "$OUT/toys.log" 2>&1 &
TOY_PID=$!
FAILED=0
wait "$DATA_PID" || FAILED=1
wait "$MC_PID" || FAILED=1
wait "$TOY_PID" || FAILED=1
[[ "$FAILED" == 0 ]] || { echo 'Ensemble failed; inspect retained failure records.'; exit 1; }
"${PYTHON_CMD}" scripts/report_preliminary_pi0.py --output "$OUT" > "$OUT/report.log" 2>&1
