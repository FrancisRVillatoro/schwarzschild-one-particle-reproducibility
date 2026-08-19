#!/usr/bin/env python3
"""Reproduce the channel-C endpoint-Hankel audit in units 2M=1.

Historical channel C: k=2, l=2, mu=2, first-order WKB preparation at
r0=1-1e-3 including the WKB amplitude-derivative term.

The terminal Bogoliubov coefficient is obtained by propagating the exact
endpoint-Hankel mode outward and projecting the WKB-prepared mode onto it
with the conserved Wronskian; no finite-radius fit of n(r) is used.
"""
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp, quad
from scipy.optimize import root
from scipy.special import hankel2
from pathlib import Path


def ffun(r): return (1-r)/r

def qfun(r): return r**1.5*np.sqrt(1-r)

def Qpump(r): return -(0.5 + 1/(16*(1-r)))/r**3

def omega2(r,k,l,mu): return k*k/ffun(r) + l*(l+1)/r**2 + mu*mu

def Omeff2(r,k,l,mu): return omega2(r,k,l,mu)-Qpump(r)

def dOmeff_dr(r,k,l,mu):
    L2=l*(l+1.0)
    domega2=k*k/(1-r)**2-2*L2/r**3
    dQ=1.5*r**(-4)-(1/16)*(-3*r**(-4)/(1-r)+r**(-3)/(1-r)**2)
    return (domega2-dQ)/(2*np.sqrt(Omeff2(r,k,l,mu)))

def ax(r): return np.sqrt(ffun(r))/r

def wkb_initial_chi(r0,k,l,mu,include_amplitude_derivative=True):
    Om=np.sqrt(Omeff2(r0,k,l,mu))
    chi=1/np.sqrt(2*Om)
    if include_amplitude_derivative:
        Omdot=-np.sqrt(ffun(r0))*dOmeff_dr(r0,k,l,mu)
        chidot=(-1j*Om-Omdot/(2*Om))*chi
    else:
        chidot=-1j*Om*chi
    return chi, -chidot/ax(r0)

def b_geometry(r,k,l):
    kap2=k*k+l*(l+1.0)
    h2=l*(l+1.0)/kap2 if kap2>0 else 0.0
    D=(1-h2)*r**4+h2*r*(1-r)
    Dp=4*(1-h2)*r**3+h2*(1-2*r)
    Dpp=12*(1-h2)*r**2-2*h2
    b=np.sqrt(D)
    L=Dp/(2*D)
    Lp=Dpp/(2*D)-Dp*Dp/(2*D*D)
    w=r*(1-r); wp=1-2*r
    s=w*(w*L*L-2*wp*L-2*w*Lp)/(4*D)
    return b,L,Lp,w,wp,s,r**3*(1-r)

def eta_tilde(r,k,l):
    if l==0:
        return -r-np.log1p(-r)
    return quad(lambda rr:b_geometry(rr,k,l)[0]/(rr*(1-rr)),0,r,
                epsabs=2e-13,epsrel=2e-13,limit=300)[0]

def integrate_endpoint_mode(k,l,mu,r_eps,r0,rtol=2e-10):
    kap=np.sqrt(k*k+l*(l+1.0)); eta=eta_tilde(r_eps,k,l); z=kap*eta
    c=np.sqrt(np.pi/4); H0=hankel2(0,z); H1=hankel2(1,z)
    v=c*np.sqrt(eta)*H0
    veta=c*(H0/(2*np.sqrt(eta))-kap*np.sqrt(eta)*H1)
    b,L,Lp,w,wp,s,q2=b_geometry(r_eps,k,l)
    vr=veta*b/w; vx=r_eps*vr
    def rhs(x,Z):
        r=np.exp(x); u=Z[0]+1j*Z[1]; ux=Z[2]+1j*Z[3]
        b,L,Lp,w,wp,s,q2=b_geometry(r,k,l); g=b/w
        Ar=L-wp/w
        F=k*k+l*(l+1.0)+s+mu*mu*q2/(b*b)
        uxx=-(-1-r*Ar)*ux-r*r*g*g*F*u
        return [ux.real,ux.imag,uxx.real,uxx.imag]
    return solve_ivp(rhs,(np.log(r_eps),np.log(r0)),
                     [v.real,v.imag,vx.real,vx.imag],method='DOP853',
                     rtol=rtol,atol=rtol*1e-2,max_step=0.10)

def historical_state_u_at_r0(k,l,mu,r0,amp_derivative=True):
    chi,chix=wkb_initial_chi(r0,k,l,mu,amp_derivative)
    b,L,Lp,w,wp,s,q2=b_geometry(r0,k,l); S=np.sqrt(b/qfun(r0))
    qlog=1.5/r0-1/(2*(1-r0)); dlogS=0.5*(L-qlog)
    u=S*chi; ur=S*(chix/r0+dlogS*chi)
    return u,ur

def exact_endpoint_projection(k,l,mu,r0=0.999,r_eps=1e-10,
                              amp_derivative=True,rtol=2e-10):
    sol=integrate_endpoint_mode(k,l,mu,r_eps,r0,rtol)
    Z=sol.y[:,-1]; v=Z[0]+1j*Z[1]; vr=(Z[2]+1j*Z[3])/r0
    u,ur=historical_state_u_at_r0(k,l,mu,r0,amp_derivative)
    g=b_geometry(r0,k,l)[0]/b_geometry(r0,k,l)[3]
    beta=(u*np.conj(vr)-ur*np.conj(v))/g/1j
    Wvv=(v*np.conj(vr)-vr*np.conj(v))/g
    return beta,Wvv

out=Path('.')
channels=[('A',2.,0,2.,.13054),('B',2.,1,2.,.00323),
          ('C',2.,2,2.,.00002),('D',2.,0,0.,.15362)]
rows=[]
for name,k,l,mu,old in channels:
    beta,W=exact_endpoint_projection(k,l,mu)
    rows.append([name,k,l,mu,beta.real,beta.imag,abs(beta)**2,old,W.imag])
repro=pd.DataFrame(rows,columns=['channel','k','l','mu','Re beta_H','Im beta_H',
    'exact |beta_H|^2','v4 rounded n_inf','Im W[v,v*]'])
repro.to_csv(out/'channel_C_reproduction_all_channels.csv',index=False)

conv=[]
for eps in [1e-6,1e-8,1e-10,1e-12]:
    beta,W=exact_endpoint_projection(2,2,2,r_eps=eps,rtol=1e-10)
    conv.append([eps,beta.real,beta.imag,abs(beta)**2,W.imag])
pd.DataFrame(conv,columns=['r_eps','Re beta_H','Im beta_H','|beta_H|^2','Im W[v,v*]']).to_csv(
    out/'channel_C_endpoint_convergence.csv',index=False)

sens=[]
for r0 in [.99,.999,.9999]:
    for amp in [False,True]:
        beta,_=exact_endpoint_projection(2,2,2,r0=r0,amp_derivative=amp)
        sens.append([r0,'WKB+amplitude derivative' if amp else 'leading WKB only',
                     beta.real,beta.imag,abs(beta)**2])
pd.DataFrame(sens,columns=['r0','preparation','Re beta_H','Im beta_H','|beta_H|^2']).to_csv(
    out/'channel_C_preparation_sensitivity.csv',index=False)

def Froot(X,r0):
    b,_=exact_endpoint_projection(X[0],2,X[1],r0=r0,r_eps=3e-11,rtol=3e-10)
    return [b.real,b.imag]
rr=root(lambda X:Froot(X,.999),[2.5057,.1556],tol=1e-10)
kstar,mustar=rr.x
bstar,_=exact_endpoint_projection(kstar,2,mustar,r0=.999,r_eps=3e-11)
h=1e-4
def B(K,M): return exact_endpoint_projection(K,2,M,r0=.999,r_eps=3e-11)[0]
dk=(B(kstar+h,mustar)-B(kstar-h,mustar))/(2*h)
dm=(B(kstar,mustar+h)-B(kstar,mustar-h))/(2*h)
det=np.linalg.det([[dk.real,dm.real],[dk.imag,dm.imag]])
pd.DataFrame([[kstar,mustar,bstar.real,bstar.imag,abs(bstar)**2,det]],columns=[
    'k*','mu*','Re beta_H','Im beta_H','|beta_H|^2','det d(Re beta,Im beta)/d(k,mu)']).to_csv(
    out/'channel_C_zero_resonance.csv',index=False)

shift=[]; x0=np.array([kstar,mustar])
for r0 in [.99,.999,.9999]:
    rr0=root(lambda X:Froot(X,r0),x0,tol=1e-9)
    b,_=exact_endpoint_projection(rr0.x[0],2,rr0.x[1],r0=r0,r_eps=3e-11)
    shift.append([r0,rr0.success,rr0.x[0],rr0.x[1],abs(b)**2]); x0=rr0.x
pd.DataFrame(shift,columns=['r0','root_success','k*','mu*','|beta_H|^2']).to_csv(
    out/'channel_C_zero_shift_with_r0.csv',index=False)
print(repro.to_string(index=False))
