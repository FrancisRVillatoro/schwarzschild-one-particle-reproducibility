"""Numerical check of Eq. (6.6) [modulus; the overall phase of beta_H relative to the Jost convention used here contains a factor exp(i kappa c(lhat)) from the constant offset between the s-wave tortoise time and the CMTV conformal time at the horizon, so only |beta_H| is compared; kappa = 80 is omitted because |beta_H|^2 ~ 1e-10 reaches the normalization error] of paper II along rays of fixed direction lhat:
   beta_H = -i/(6 pi lhat^2 kappa^2) + o(kappa^-2),  |beta_H|^2 = 1/(36 pi^2 lhat^4 kappa^4).
Exact horizon-Jost solution with K = kappa sqrt(1-lhat^2), L^2 = kappa^2 lhat^2 (real L), projected on the endpoint pair."""
import numpy as np
from num_checks import propagate, EG
def coeffs(K,L2,xmin=np.log(1e-12)):
    ph,px=propagate(K,L2,np.exp(1j*np.pi/4),0j,xmin=xmin)
    b=px; a=ph-b*xmin; L=np.sqrt(L2)
    av=np.sqrt(np.pi/2)*(1-2j/np.pi*(np.log(L)+EG)); bv=-1j/np.sqrt(2*np.pi)
    return np.linalg.solve(np.array([[av,np.conj(av)],[bv,np.conj(bv)]]),np.array([a,b]))
for lh2 in (0.5,0.8):
    print(f'\nlhat^2 = {lh2}:   kappa   beta_H (numerical)        -i/(6 pi lhat^2 kappa^2)    |beta|^2 kappa^4 (num)   1/(36 pi^2 lhat^4)')
    for kap in (10.,20.,40.):
        K=kap*np.sqrt(1-lh2); L2=kap**2*lh2
        al,be=coeffs(K,L2)
        print(f'               {kap:5.0f}   {be.real:+.3e}{be.imag:+.3e}i   {-1/(6*np.pi*lh2*kap**2):+.3e}i   {abs(be)**2*kap**4:.5f}   {1/(36*np.pi**2*lh2**2):.5f}   norm-1={abs(al)**2-abs(be)**2-1:+.0e}')
