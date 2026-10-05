#!/usr/bin/env python3
"""
Check that the shrinking infrared layer lambda=O(mu) of the exact massive
Kasner modes leaves no finite contribution to <Phi^2>_ren as mu -> 0.

The only finite residue comes from the second-order adiabatic counterterm
    -mu^4/[16 (lambda^2+mu^2)^(7/2)],
whose subtraction produces +1/(240 pi^2).

Requires numpy, scipy, pandas.
"""
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from scipy.special import hankel2
from numpy.polynomial.legendre import leggauss

def Gd(s,L,u,mu):
    a=L*L*u*u*np.exp(8*s/3)
    b=L*L*(1-u*u)*np.exp(2*s/3)
    c=mu*mu*np.exp(2*s)
    G=a+b+c
    G1=(8*a+2*b)/3+2*c
    G2=(64*a+4*b)/9+4*c
    G3=(512*a+8*b)/27+8*c
    return G,G1,G2,G3

def W2(s,L,u,mu):
    G,G1,G2,G3=Gd(s,L,u,mu)
    H=G-G2/(4*G)+5*G1*G1/(16*G*G)
    H1=G1-.25*(G3/G-G2*G1/G**2)+(5/8)*(G1*G2/G**2-G1**3/G**3)
    W=np.sqrt(H)
    return W,H1/(2*W)

def start(L,u,mu,target=70.):
    def Om(s): return np.sqrt(Gd(s,L,u,mu)[0])
    if Om(0)>=target: return 0.
    hi=1.
    while Om(hi)<target: hi*=1.5
    return brentq(lambda s:Om(s)-target,0,hi)

def amp2(L,u,mu,target=70.,rtol=2e-9):
    S=start(L,u,mu,target)
    W,W1=W2(S,L,u,mu)
    R0=-W1/(2*W)-1j*W
    if S==0: return -1/(2*R0.imag)
    def rhs(s,z):
        R=z[0]+1j*z[1]
        d=-Gd(s,L,u,mu)[0]-R*R
        return [d.real,d.imag]
    sol=solve_ivp(rhs,(S,0),[R0.real,R0.imag],method="DOP853",
                  rtol=rtol,atol=rtol*1e-2)
    R=sol.y[0,-1]+1j*sol.y[1,-1]
    return -1/(2*R.imag)

def ir_exact_difference(mu,X=8.,Nx=14,Nu=14):
    xx,wx=leggauss(Nx); xs=X*(xx+1)/2; wx=X*wx/2
    uu,wu=leggauss(Nu); us=(uu+1)/2; wu=wu/2
    tot=0.
    for x,wxj in zip(xs,wx):
        L=mu*x
        av=0.
        for u,wuj in zip(us,wu):
            av += wuj*(amp2(L,float(u),mu)-amp2(L,float(u),0.0))
        tot += wxj*x*x*av
    return mu**3*tot/(2*np.pi**2)

def counterterm_parts(mu,X=8.):
    ordinary=quad(
        lambda x:x*x*(-(1/(2*np.sqrt(1+x*x))-1/(2*x))),
        0,X,epsabs=1e-12,epsrel=1e-12
    )[0]*mu**2/(2*np.pi**2)
    residue=quad(
        lambda x:x*x/(1+x*x)**3.5,
        0,X,epsabs=1e-12,epsrel=1e-12
    )[0]/(32*np.pi**2)
    return ordinary,residue

# Exact lambda=0 check
for mu in [0.2,0.05,0.01]:
    num=amp2(0,0.3,mu,target=80,rtol=1e-10)
    exact=(np.pi/4)*abs(hankel2(0,mu))**2
    print("axis",mu,num,exact,num/exact-1)

rows=[]
for mu in [0.1,0.05,0.025,0.0125,0.00625]:
    D=ir_exact_difference(mu)
    O,R=counterterm_parts(mu)
    rows.append([mu,D,D/(mu**3*(1+np.log(mu)**2)),O,R,D+O+R])
df=pd.DataFrame(rows,columns=[
    "mu","exact_mode_difference","scaled_exact_difference",
    "ordinary_counterterm","residue_X8","renormalized_layer_difference"])
print(df.to_string(index=False))
print("full residue =",1/(240*np.pi**2))
