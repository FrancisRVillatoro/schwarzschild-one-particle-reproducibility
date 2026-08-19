#!/usr/bin/env python3
"""
Reproduce the Kasner vacuum-polarization coefficient used in the manuscript audit.

Requires: numpy, scipy

Mode equation in s = ln y:
    psi_ss + lambda^2 [u^2 exp(8s/3)+(1-u^2)exp(2s/3)] psi = 0

The positive-frequency solution is imposed at large s by second-order WKB.
The Riccati variable R=psi_s/psi obeys
    R_s = -Omega^2 - R^2,
and Wronskian normalization gives
    |psi(0)|^2 = -1/(2 Im R(0)).

After angular integration the massless second-order adiabatic counterterm is
exactly 1/(2 lambda), because its lambda^-3 angular coefficient integrates
to zero.  The nonuniform massless-limit residue is 1/(240 pi^2).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss

def Fd(s,u):
    a=u*u*np.exp(8*s/3); b=(1-u*u)*np.exp(2*s/3)
    return a+b,(8*a+2*b)/3,(64*a+4*b)/9,(512*a+8*b)/27

def W2d(s,L,u):
    F,F1,F2,F3=Fd(s,u)
    G=L*L*F-F2/(4*F)+5*F1*F1/(16*F*F)
    G1=L*L*F1-.25*(F3/F-F2*F1/F**2)+(5/8)*(F1*F2/F**2-F1**3/F**3)
    W=np.sqrt(G)
    return W,G1/(2*W)

def amp2(L,u,target=60.,rtol=3e-9):
    targ=max(target,2*L)
    def h(s): return L*np.sqrt(Fd(s,u)[0])-targ
    if h(0)>=0: S=0.
    else:
        hi=1.
        while h(hi)<0: hi*=2
        S=brentq(h,0,hi)
    W,W1=W2d(S,L,u)
    R0=-W1/(2*W)-1j*W
    if S==0: return -1/(2*R0.imag)
    def rhs(s,z):
        R=z[0]+1j*z[1]
        d=-L*L*Fd(s,u)[0]-R*R
        return [d.real,d.imag]
    sol=solve_ivp(rhs,(S,0),[R0.real,R0.imag],
                  method="DOP853",rtol=rtol,atol=rtol*1e-2)
    R=sol.y[0,-1]+1j*sol.y[1,-1]
    return -1/(2*R.imag)

def g(L,Nu=24,target=60.,rtol=3e-9):
    x,w=leggauss(Nu); us=(x+1)/2; ws=w/2
    A=sum(ww*amp2(L,float(u),target,rtol) for u,ww in zip(us,ws))
    return L*L*(A-1/(2*L))

def segment(a,b,n=16,Nu=24,target=60.,rtol=3e-9):
    x,w=leggauss(n); L=(a+b)/2+(b-a)*x/2
    return (b-a)/2*np.dot(w,[g(float(z),Nu,target,rtol) for z in L])

raw=sum(segment(a,b) for a,b in [(0,.25),(.25,1),(1,3),(3,8)])
C8=raw/(2*np.pi**2)

# Sixth-order adiabatic expansion fixes the first omitted angular coefficient
# exactly:
#   int_0^1 |psi|^2 du - 1/(2L) = c7 L^-7 + O(L^-9),
#   c7 = -1/2835.
# Since g(L)=L^2[...], its leading term is c7 L^-5.  We fit only the
# subleading L^-7 and L^-9 terms on a numerically safe window.
c7_exact=-1/2835
Ls=np.array([6.,8.,10.,12.])
gs=np.array([g(L,Nu=28,target=250.,rtol=3e-10) for L in Ls])
res=gs-c7_exact*Ls**-5
M=np.column_stack([Ls**-7,Ls**-9])
c9,c11=np.linalg.lstsq(M,res,rcond=None)[0]
tail=(c7_exact/(4*8**4)+c9/(6*8**6)+c11/(8*8**8))/(2*np.pi**2)

Cdyn=C8+tail
CIR=1/(240*np.pi**2)
CK=Cdyn+CIR

qJ=7.1782e-4
CJ=2*qJ/9

print("C_dyn[0,8] =",C8)
print("analytic c7   =",c7_exact)
print("tail          =",tail)
print("C_dyn total   =",Cdyn)
print("1/(240 pi^2)  =",CIR)
print("C_K           =",CK)
print("C_J           =",CJ)
print("relative diff =", (CK-CJ)/CJ)
print("q_K           =",4.5*CK)
