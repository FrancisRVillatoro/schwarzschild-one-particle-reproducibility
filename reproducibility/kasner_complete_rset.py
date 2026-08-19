#!/usr/bin/env python3
"""
Reproduce the leading terminal Kasner RSET for the massless minimally coupled
scalar scaling state used in fate_manuscript_v7.1.

Requires kasner_rset_core.py in the same directory.

Main idea:
1. Compute the single independent anisotropic coefficient
      Delta = C_parallel - C_perp
   from the exact massless Kasner modes.
2. Use the angularly integrated fourth-order adiabatic counterterm
      J_0(lambda) = int_0^1 (3u^2-1) A_ad^(4) du = -1/(63 lambda^3).
3. Add the nonuniform small-mass residue
      Delta_IR = -523/(498960 pi^2).
4. Reconstruct all diagonal RSET coefficients from exact conservation
   and the trace identity using the independently determined
      C_K = 1.59555e-4
   for <Phi^2>_ren = C_K sigma^-2.
"""
import numpy as np
from pathlib import Path
from kasner_rset_core import integrate_delta_massless_interval

CK = 1.59555e-4
intervals = [(0.0,0.2),(0.2,1.0),(1.0,4.0),(4.0,16.0)]

pieces=[]
for lo,hi in intervals:
    val=integrate_delta_massless_interval(lo,hi,Nu=26,Nl=18)
    pieces.append(val)
    print(f"{lo:8g} {hi:8g} {val:.12e}")

# Sixth-order adiabatic expansion fixes the first omitted coefficient
# exactly:
#   B_Delta(L) = int(3u^2-1)|psi|^2 du + 1/(63 L^3)
#              = d7 L^-7 + d9 L^-9 + ... ,
#   d7 = 4/2835.
# The radial integrand is L^4 B_Delta/(4 pi^2).  We use the exact d7
# tail and estimate d9 from the numerically safe [4,16] panel. Higher
# terms are covered by the quoted 3e-9 absolute uncertainty in Delta_dyn.
d7=4/2835
a,b=4.0,16.0
lead_panel=d7/(8*np.pi**2)*(a**-2-b**-2)
coef_d9=(a**-4-b**-4)/(16*np.pi**2)
d9_est=(pieces[-1]-lead_panel)/coef_d9
tail_d7=d7/(8*np.pi**2*b**2)
tail_d9=d9_est/(16*np.pi**2*b**4)
tail=tail_d7+tail_d9
Delta_dyn=sum(pieces)+tail
Delta_dyn_unc=3e-9

# Analytic nonuniform fourth-order residue:
Delta_IR=-523/(498960*np.pi**2)
Delta=Delta_dyn+Delta_IR

A_geom=1/(1215*np.pi**2)
Trace=2*CK+A_geom

C_rho=(Trace-2*Delta)/8
C_perp=3*C_rho+Delta/3
C_parallel=C_perp+Delta

lam_t=(81/4)*C_parallel
lam_r=-(81/4)*C_rho
lam_perp=(81/4)*C_perp
qK=(9/2)*CK

print("\nd7 exact   =",d7)
print("d9 estimate=",d9_est)
print("UV tail    =",tail)
print("Delta_dyn  =",Delta_dyn,"+/-",Delta_dyn_unc)
print("Delta_IR  =",Delta_IR)
print("Delta_K   =",Delta)
print("\nKasner coefficients:")
print("C_rho      =",C_rho)
print("C_parallel =",C_parallel)
print("C_perp     =",C_perp)
print("\nSchwarzschild coefficients:")
print("lambda_t    =",lam_t)
print("lambda_r    =",lam_r)
print("lambda_perp =",lam_perp)
print("q_K         =",qK)
