"""Phase 2 numerics for the fermionic paper III (units 2M = 1; modes started in the instantaneous positive-energy
eigenspace at r0 = 0.9; conventions of the manuscript).  Theorem B: convergence of the Parker count n_P and of the
reference count n_ref to the terminal occupation n_perp, with the local exponents of |n - n_perp| between successive radii."""
import numpy as np
from dirac_phase0 import H, proj, g0, g2, evolve
rs=np.array([1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-9])
for (k,kap,m) in [(5,2,0.),(30,3,0.),(100,1,0.),(5,2,1.)]:
    Ys=evolve(k,kap,m,0.9,rs); Pm=proj(g0@g2*np.sign(kap),-1)
    nperp=np.trace(Ys[-1].conj().T@Pm@Ys[-1]).real/2
    nP=[np.trace(Y.conj().T@proj(H(r,k,kap,m),-1)@Y).real/2 for r,Y in zip(rs,Ys)]; nR=[np.trace(Y.conj().T@Pm@Y).real/2 for Y in Ys]
    print(f'mode (k,lambda,m)=({k},{kap},{m}):  n_perp = {nperp:.8f}')
    print('      r         n_P           n_ref        |n_P-n_perp|   exp    |n_ref-n_perp|   exp')
    for i,r in enumerate(rs[:-1]):
        dP,dR=abs(nP[i]-nperp),abs(nR[i]-nperp)
        eP=np.log(abs(nP[i]-nperp)/max(abs(nP[i+1]-nperp),1e-16))/np.log(rs[i]/rs[i+1]) if i<len(rs)-2 else float('nan')
        eR=np.log(abs(nR[i]-nperp)/max(abs(nR[i+1]-nperp),1e-16))/np.log(rs[i]/rs[i+1]) if i<len(rs)-2 else float('nan')
        print(f'  {r:8.0e}   {nP[i]:.8f}   {nR[i]:.8f}   {dP:.2e}   {eP:5.2f}   {dR:.2e}     {eR:5.2f}')
