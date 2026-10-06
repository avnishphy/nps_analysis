#!/usr/bin/env bash
# Source the Hall C/NPS environment before invoking.
set -euo pipefail
src="$(cd "$(dirname "$0")/.." && pwd)"
config="${1:?usage: run_sigparam_validation.sh CONFIG NEW_OUT}"
out="${2:?provide a new output directory}"
if [[ -e "$out" ]]; then echo "Output exists: $out" >&2; exit 1; fi
mkdir -p "$out/build";out="$(cd "$out" && pwd)"
python3 "$src/generate_xsec_config.py" "$config" "$out/build/xsec_config.h"
read -r -a flags <<< "$(root-config --cflags --libs)"
for test in test_sigparam_events validate_vertex_epsilon test_proxy_staged; do
  g++ "$src/tests/$test.C" -I"$src" -I"$out/build" "${flags[@]}" -lMinuit2 -O2 -std=c++17 -o "$out/build/$test"
done
"$out/build/test_sigparam_events" | tee "$out/event_test.log"
"$out/build/test_proxy_staged" | tee "$out/staged_cone_test.log"
python3 -m unittest discover -s "$src/tests" -p 'test_*.py'
for method in simc_model no_simc_model; do
  g++ "$src/excl_xsec_pi0_analysis_${method}.C" -I"$src" -I"$out/build" "${flags[@]}" -lMinuit2 -O2 -std=c++17 -o "$out/build/$method"
done
g++ -O2 -std=c++17 -I"$src" "$src/tests/sigparam_evaluate.cpp" -o "$out/build/sigparam_evaluate"
bash "$src/tests/run_sigparam_fortran.sh" /u/group/nps/singhav/simc_gfortran_updated/physics_pion.f "$out/fortran"
for kind in nominal matching-failures empty-row; do
  variant=();if [[ "$kind" != nominal ]]; then variant=("$kind"); fi
  python3 "$src/tests/make_sigparam_fixture.py" "$config" "$out/$kind" "$out/build/sigparam_evaluate" "${variant[@]}"
  objectives=(gaussian);if [[ "$kind" != matching-failures ]]; then objectives+=(scaled-poisson); fi
  for objective in "${objectives[@]}"; do
    "$out/build/simc_model" --data-file "$out/$kind/data.root" --sim-file "$out/$kind/sim.root" --vertex_simc_file "$out/$kind/raw.root" --out-dir "$out/${kind}_${objective}" --fit-objective "$objective" --no-png --no-pdf --no-diagnostics > "$out/${kind}_${objective}.log" 2>&1
  done
done
python3 - "$out" <<'PY'
import csv,json,sys
from pathlib import Path
import numpy as np
p=Path(sys.argv[1]);truth=np.array(json.loads((p/'nominal/truth.json').read_text())['theta'])
for mode in ('gaussian','scaled-poisson'):
    rows=list(csv.DictReader((p/f'nominal_{mode}/model_parameters.csv').open()))
    assert [r['name'] for r in rows]==['N_U','DeltaB_U','N_LT','N_TT']
    fitted=np.array([float(r['value']) for r in rows]);error=max(abs(fitted-truth))
    assert error<1e-4,error
    print('CLOSURE',mode,'max_parameter_absolute',error)
r=next(csv.DictReader((p/'matching-failures_gaussian/model_vertex_matching.csv').open()))
assert all(int(r[k])==1 for k in ('unmatched','duplicate','sigcm_mismatch','invalid_vertex')),r
r=next(csv.DictReader((p/'empty-row_scaled-poisson/scaled_poisson_rows.csv').open()))
assert r['included']=='1' and r['zero_bin']=='1' and float(r['mu'])>0,r
print('PASS: matching guards and supported empty-row Poisson treatment')
PY
for objective in gaussian scaled-poisson; do
  python3 "$src/tests/validate_sigparam_exports.py" "$out/nominal_$objective"
done
echo "Validation passed: $out"
