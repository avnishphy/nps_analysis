#!/usr/bin/env bash
# Source Hall C setup before invoking. Read the actual unchanged Fortran source.
set -euo pipefail
src="$(cd "$(dirname "$0")/.." && pwd)"
original="${1:-/u/group/nps/singhav/simc_gfortran_updated/physics_pion.f}"
out="${2:?usage: bash run_sigparam_fortran.sh ORIGINAL_PHYSICS_PION_F NEW_BUILD_DIR}"
mkdir -p "$out"
python3 - "$original" "$out/sigparam_original.f" <<'PY'
from pathlib import Path
import sys
s=Path(sys.argv[1]).read_text()
start=s.index('      real*8 function sig_param_2021(')
# Both routines are the final routines in the audited source. Preserve text.
Path(sys.argv[2]).write_text(s[start:])
PY
sha256sum "$original" > "$out/fortran_source.sha256"
gfortran -O -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8 -c "$out/sigparam_original.f" -o "$out/sigparam_original.o"
g++ -std=c++17 -O2 -I"$src" "$src/tests/test_sigparam_fortran.cpp" "$out/sigparam_original.o" -lgfortran -lquadmath -o "$out/test_sigparam_fortran"
"$out/test_sigparam_fortran" | tee "$out/charged_reproduction.log"
