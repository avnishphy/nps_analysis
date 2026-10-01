#!/usr/bin/env bash
# Prerequisite: source Hall C/NPS ROOT environment through csh.
# Compare to the original SIMC Fortran and run a synthetic extraction closure.
set -euo pipefail
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
out="${1:-$(mktemp -d /tmp/nps_simc_phi_average.XXXXXX)}"
simc_source="${2:-/u/group/nps/singhav/simc_gfortran_updated/physics_pion.f}"
mkdir -p "$out"
read -r -a root_flags <<< "$(root-config --cflags --libs)"
g++ -O2 -std=c++17 "$repo_root/tests/test_simc_phi_average.C" "${root_flags[@]}" -o "$out/test_model"
"$out/test_model" "$out"
python3 - "$simc_source" "$out" <<'PY'
from pathlib import Path
import sys
source, out = Path(sys.argv[1]), Path(sys.argv[2])
text = source.read_text()
(out/'reference.f').write_text(text[text.index('      real*8 function sig_param_2021'):])
(out/'driver.f90').write_text('''program oracle
implicit none
real(kind=8) :: th,phi,t,q2,s,eps,sig_param_2021
integer :: ios
do
 read(*,*,iostat=ios) th,phi,t,q2,s,eps
 if(ios /= 0) exit
 write(*,'(ES26.17)') sig_param_2021(th,phi,t,q2,s,eps,0,.true.)
end do
end program
''')
PY
gfortran -O2 -fdefault-real-8 -ffixed-line-length-132 "$out/reference.f" "$out/driver.f90" -o "$out/oracle"
"$out/oracle" < "$out/model_cases.txt" > "$out/model_fortran.txt"
extractor="$(python3 "$repo_root/tests/build_xsec_test_binary.py" \
    --out-dir "$out/build" --phi-bins 12 --t-edges=-2,0 \
    --q2-edges=5.1,6.4 --xb-edges=.45,.65)"
"$extractor" --data-file "$out/data.root" --sim-file "$out/sim.root" \
    --out-dir "$out/extracted" --target-contam 1 --target-contam-err 0 --fixed-default-model \
    --no-png --no-pdf --no-diagnostics > "$out/extraction.log" 2>&1
"$extractor" --data-file "$out/data.root" --sim-file "$out/sim.root" \
    --out-dir "$out/scaled" --target-contam 1 --target-contam-err 0 --fixed-default-model \
    --simc-yield-scale 2 --no-png --no-pdf --no-diagnostics > "$out/scaled.log" 2>&1
python3 - "$out" <<'PY'
import csv, math, sys
from pathlib import Path
out=Path(sys.argv[1])
cpp=[float(x) for x in (out/'model_cpp.txt').read_text().split()]
fort=[float(x) for x in (out/'model_fortran.txt').read_text().split()]
assert len(cpp)==len(fort)==324
for a,b in zip(cpp,fort):
    assert math.isclose(a,b,rel_tol=2e-12,abs_tol=1e-22), (a,b)
rows=list(csv.DictReader((out/'extracted/excl_xsec_pi0_analysis_simc_model_summary.csv').open()))
assert len(rows)==12
for row in rows:
    assert math.isclose(float(row['ratio']),1.7,rel_tol=5e-6)
    assert math.isclose(float(row['xsec']),float(row['ratio'])*float(row['model_xsec_phi_center']),rel_tol=1e-8)
for key in ('q2_ref','xb_ref','tprime_ref','t_ref','W_ref'):
    assert len({r[key] for r in rows})==1, key
assert any(not math.isclose(float(r['model_sigcm_full_weight_mean']),
                           float(r['model_xsec_phi_center']),rel_tol=1e-3) for r in rows)
slices=list(csv.DictReader((out/'extracted/excl_xsec_pi0_analysis_simc_model_slice_summary.csv').open()))
assert slices[0]['fit_xsec_ok']=='1'
assert float(slices[0]['fit_xsec_chi2'])<1e-11
assert slices[0]['sigmaTLp_available']=='0'
assert slices[0]['sigmaTLp_status']=='unavailable_missing_helicity_luminosities_and_beam_polarization'
assert math.isnan(float(slices[0]['sigmaTLp']))
assert math.isnan(float(slices[0]['sigmaTLp_err']))
scaled=list(csv.DictReader((out/'scaled/excl_xsec_pi0_analysis_simc_model_summary.csv').open()))
for a,b in zip(rows,scaled):
    assert math.isclose(float(a['model_xsec_phi_center']),float(b['model_xsec_phi_center']),rel_tol=1e-9)
    assert math.isclose(float(a['xsec']),2*float(b['xsec']),rel_tol=1e-9)
    assert math.isclose(float(a['xsec_err']),2*float(b['xsec_err']),rel_tol=1e-9)
print('PASS 324 original-Fortran comparisons and 12-bin extraction/fit closure')
print('PASS physical SIMC normalization scaling')
print("PASS SIMC-model LT' unavailable output contract")
print(f'Artifacts: {out}')
PY
