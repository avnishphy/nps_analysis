#!/usr/bin/env python3
"""Compare the default C++ pi0 evaluator to the local authoritative SIMC Fortran."""
import math
from pathlib import Path
import subprocess
import sys
import tempfile

repo = Path(__file__).resolve().parents[1]
source = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/u/group/nps/singhav/simc_gfortran_updated/physics_pion.f")
text = source.read_text()
tail = text[text.index("      real*8 function sig_param_2021"):]
points = [(3.7, 2.813, -.2218935, .1, .73),
          (4.5, 2.75, -.5, 1.3, .7),
          (3.3, 2.9, -.35, 4.2, .75),
          (4.7, 2.7, -.8, 3.14, .71)]
mp, mpi = .9382720813, .1349768
lines = []
for q2, w, t, phi, eps in points:
    ws = w*w
    eg = (ws-mp*mp-q2)/(2*w)
    pg = math.sqrt(eg*eg+q2)
    epi = (ws+mpi*mpi-mp*mp)/(2*w)
    pp = math.sqrt(epi*epi-mpi*mpi)
    theta = math.acos((t+q2-mpi*mpi+2*eg*epi)/(2*pg*pp))
    lines.append(f"{q2:.16g} {w:.16g} {t:.16g} {theta:.16g} {phi:.16g} {eps:.16g}")
cases = "\n".join(lines)+"\n"
with tempfile.TemporaryDirectory(prefix="nps_fortran_model_") as td:
    tmp = Path(td)
    (tmp/"oracle.f").write_text(tail + """
      program oracle
      implicit none
      real*8 q2,w,t,th,phi,eps,v,sig_param_2021
      integer ios
 10   read(*,*,iostat=ios) q2,w,t,th,phi,eps
      if (ios.ne.0) goto 20
      v=sig_param_2021(th,phi,t,q2,w*w,eps,0,.true.)
      write(*,'(ES24.16)') v
      goto 10
 20   continue
      end
""")
    (tmp/"model.cpp").write_text("""
#include <iostream>
#include "src/xsec_extract/simc_pi0_reweight.h"
int main() {
 double q2,w,t,theta,phi,eps;
 while (std::cin>>q2>>w>>t>>theta>>phi>>eps) {
  std::cout.precision(17);
  std::cout<<nps_pi0_reweight::evaluate(q2,w,t,phi,eps,
   .9382720813,.1349768,nps_pi0_reweight::defaults())<<"\\n";
 }
}
""")
    subprocess.run(["gfortran", "-O2", "-fdefault-real-8", "-ffixed-line-length-132",
                    str(tmp/"oracle.f"), "-o", str(tmp/"oracle")], check=True)
    subprocess.run(["g++", "-O2", "-std=c++17", "-I"+str(repo),
                    str(tmp/"model.cpp"), "-o", str(tmp/"model")], check=True)
    a = [float(x) for x in subprocess.check_output([str(tmp/"oracle")],
                                                     input=cases.encode()).split()]
    b = [float(x) for x in subprocess.check_output([str(tmp/"model")],
                                                     input=cases.encode()).split()]
    assert len(a) == len(b) == len(points)
    max_rel = max(abs(x/y-1) for x, y in zip(a, b))
    assert max_rel < 1e-12, (a, b, max_rel)
    print(f"PASS {len(points)} Fortran/C++ physical points; max relative difference {max_rel:.3g}")
