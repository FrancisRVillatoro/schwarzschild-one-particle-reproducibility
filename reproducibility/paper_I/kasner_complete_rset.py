#!/usr/bin/env python3
"""Final Kasner coefficients and every number of Table 1 and Section 5 of paper I.

Delta_dyn: baseline mode integrals on the panels (0,0.2,1,4,16) with Gauss-Legendre
quadrature (Nu = 26 angular nodes, Nl = 18 radial nodes per panel, Liouville-Green data
of fourth adiabatic order at local frequency 120, rtol 1e-10), plus the ultraviolet tail
beyond lambda = 16 obtained from a direct fit of the bracket
    B(lambda) = int_0^1 (3u^2-1)|psi|^2 du + 1/(63 lambda^3) = d7/lambda^7 + d9/lambda^9 + d11/lambda^11 + ...
on lambda in [8,16] with d7 = 4/2835 exact (this supersedes the earlier estimate of d9
from the last panel alone, which neglected d11). Systematic shifts measured in the
convergence study (convergence_results.txt) are applied as corrections.

C_dyn: baseline value of reproduce_kasner_phi2.py plus the systematic corrections of the
convergence study (initial time, angular quadrature with its geometric extrapolation,
radial quadrature). Requires kasner_rset_core.py and convergence_study.py in the same
directory."""
import numpy as np, sys, time
sys.argv=['x','none']
exec(open('convergence_study.py').read().split("which=sys.argv[1]")[0])   # exact_A, Dpiece
from numpy.polynomial.legendre import leggauss
t0=time.time()
# ---------------- Delta_dyn ----------------
edges=[0.0,0.2,1.0,4.0,16.0]
pieces=[Dpiece(a,b,26,18,120.,1e-10) for a,b in zip(edges[:-1],edges[1:])]
print('Delta_dyn panels:',[f'{p:.10e}' for p in pieces],f'({time.time()-t0:.0f}s)')
d7=4/2835
def Bbr(L,Nu=48,target=360.,rtol=1e-11):
    xu,wu=leggauss(Nu); us=(xu+1)/2; wu=wu/2
    return sum(wU*(3*u*u-1)*exact_A(L,float(u),target,rtol) for u,wU in zip(us,wu))+1/(63*L**3)
Ls=np.arange(8.,16.01,1.); yfit=np.array([Bbr(L)*L**7-d7 for L in Ls])
tails={}
for name,mask,ks in [('d9,d11 on [8,16]',Ls>=8,(9,11)),('d9,d11,d13 on [8,16]',Ls>=8,(9,11,13)),('d9,d11 on [10,16]',Ls>=10,(9,11))]:
    Mfit=np.column_stack([Ls[mask]**(7-k) for k in ks]); cf=np.linalg.lstsq(Mfit,yfit[mask],rcond=None)[0]
    dk={7:d7}; dk.update(dict(zip(ks,cf)))
    tails[name]=sum(v/((k-5)*16.**(k-5)) for k,v in dk.items())/(4*np.pi**2)
    print(f'tail beyond 16, fit {name}: {tails[name]:.6e}   coefficients {dict((k,round(float(v),5)) for k,v in dk.items())}')
tail=np.mean(list(tails.values())); tail_unc=(max(tails.values())-min(tails.values()))/2
dsys_D=-6.0e-11-1.1e-11          # rtol 1e-10 -> 1e-12 and initial time 120 -> 240 (convergence_results.txt)
Delta_dyn=sum(pieces)+tail+dsys_D
Delta_dyn_unc=np.sqrt(tail_unc**2+(6e-11)**2+(2e-11)**2+(3e-12)**2)   # tail spread, tolerance, cutoff check, radial nodes
Delta_IR=-523/(498960*np.pi**2); Delta_K=Delta_dyn+Delta_IR
# ---------------- C_dyn ----------------
C_dyn_base=-2.62617123305349e-4   # reproduce_kasner_phi2.py baseline (Nu=24, n=16, target 60, rtol 3e-9, cutoff 8, fit on [6,12])
corr_time=+1.23e-10; corr_ang=-5.15e-10-0.40e-10; corr_rad=+0.45e-10   # convergence_results.txt (angular: 24->40 plus geometric extrapolation)
C_dyn=C_dyn_base+corr_time+corr_ang+corr_rad
C_dyn_unc=np.sqrt((4.06e-10)**2+(0.4e-10)**2)   # cutoff 8->6 (tail test), angular extrapolation
C_IR=1/(240*np.pi**2); C_K=C_dyn+C_IR
# ---------------- derived coefficients (conservation 4C_perp = C_par + 9C_rho, trace -C_rho + C_par + 2C_perp = 2C_K + A) ----------------
A_geom=1/(1215*np.pi**2); T=2*C_K+A_geom
C_perp=3/8*(T-10/9*Delta_K); C_par=C_perp+Delta_K; C_rho=(3*C_perp-Delta_K)/9
lam_t,lam_r,lam_perp,q=81/4*C_par,-81/4*C_rho,81/4*C_perp,4.5*C_K
print(f'\nC_dyn = {C_dyn:.7e} +/- {C_dyn_unc:.1e}      C_K = {C_K:.7e}      q_K = {q:.6e}')
print(f'Delta_dyn = {Delta_dyn:.7e} +/- {Delta_dyn_unc:.1e}   Delta_K = {Delta_K:.7e}')
print(f'minimal coupling: C_rho, C_par, C_perp = {C_rho:.6e} {C_par:.6e} {C_perp:.6e}')
print(f'                  lambda_t, lambda_r, lambda_perp = {lam_t:.6e} {lam_r:.6e} {lam_perp:.6e}   Delta (4/81)(lt-lp) = {4/81*(lam_t-lam_perp):.6e}')
xi=1/6
print(f'conformal coupling: C = {C_rho-2*xi*C_K:.4e} {C_par-10/3*xi*C_K:.4e} {C_perp-16/3*xi*C_K:.4e}   lambda = {81/4*(C_par-10/3*xi*C_K):.4e} {-81/4*(C_rho-2*xi*C_K):.4e} {81/4*(C_perp-16/3*xi*C_K):.4e}')
print(f'null-projection thresholds: xi_par = {3*(C_rho+C_par)/(16*C_K):.5f}  xi_perp = {3*(C_rho+C_perp)/(22*C_K):.5f}  xi_rho = {C_rho/(2*C_K):.5f}')
print(f'trace check: -C_rho+C_par+2C_perp - (2C_K+A) = {-C_rho+C_par+2*C_perp-T:.1e};  conformal trace - A = {-(C_rho-2*xi*C_K)+(C_par-10/3*xi*C_K)+2*(C_perp-16/3*xi*C_K)-A_geom:.1e}')
qJ,ltJ,lrJ,lpJ=7.1782e-4,2.6357e-3,-1.1926e-3,3.3494e-3
print(f'Jiang-Jiang comparison: (q-qJ)/qJ = {(q-qJ)/qJ:.2e}   Delta_J = {4/81*(ltJ-lpJ):.4e}   relative difference in Delta = {(Delta_K-4/81*(ltJ-lpJ))/Delta_K:.2e}   ({time.time()-t0:.0f}s)')
