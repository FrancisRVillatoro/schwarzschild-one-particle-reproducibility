"""Independent universal transition computation with superadiabatic Jost starting data.
Does not treat floating-point convergence as a certified interval bound.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gamma
import sympy as sp

x, lam = sp.symbols('x lam', real=True, positive=True)
nu=sp.sqrt(1+x**6)
g= -3*x*x/(2*(1+x**6))
L=1/(4*sp.I*lam*nu)
c0=-L*g
c1=-L*(g+sp.diff(c0,x)+g*c0*c0)
c2=-L*(g+sp.diff(c1,x)+g*c1*c1)
# evaluate recursively via lambdified expression
funcs = [sp.lambdify((x,lam), c, modules='numpy', cse=True) for c in (c0,c1,c2)]

def transition(Lambda:float,rho0:float|None=None,rtol:float=3e-12,order:int=2):
    if rho0 is None: rho0=max(20.0,32/np.sqrt(Lambda))
    x0=np.sqrt(rho0)
    theta=np.arctan2(1.,x0**3)
    up=np.array([np.cos(theta/2),np.sin(theta/2)],complex)
    um=np.array([-np.sin(theta/2),np.cos(theta/2)],complex)
    c=complex(funcs[order](x0,Lambda)) if order>=0 else 0
    init=(up+c*um)/np.sqrt(1+abs(c)**2)
    def rhs(t,u):
        p,q=u
        return 2j*Lambda*np.array([t**3*p+q,p-t**3*q])
    sol=solve_ivp(rhs,(x0,0.),init,method='DOP853',rtol=rtol,atol=max(1e-15,rtol*0.01))
    if not sol.success: raise RuntimeError(sol.message)
    end=sol.y[:,-1]; norm=float(np.vdot(end,end).real)
    F=float(abs(end[0]-end[1])**2/(2*norm))
    return dict(Lambda=Lambda,rho0=rho0,rtol=rtol,order=order,F=F,nfev=sol.nfev,norm_error=norm-1,initial_c_abs=abs(c))

if __name__=='__main__':
    import time
    t0=__import__('time').time()
    for Lm in [.001,.01,.02,.1, 4.4,4.6,5,6,8,10,12]:
        for order in [-1,0,1,2]:
            r=transition(Lm,order=order)
            print(f'{Lm:8g} {order:2d} R={r["rho0"]:.1f} F={r["F"]:.12g} ne={r["nfev"]} c={r["initial_c_abs"]:.3g}',flush=True)
    print('total secs',__import__('time').time()-t0)
