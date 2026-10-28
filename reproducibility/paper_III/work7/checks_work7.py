"""Independent Work7 checks: Clifford symmetries, massive angle and operator-bound scaling.
These floating-point experiments are diagnostics, not interval-certified proofs.
"""
from __future__ import annotations
import json
import math
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

OUT=Path(__file__).resolve().parent
I=np.eye(4,dtype=complex)
gamma0=np.diag([1.,1.,-1.,-1.]).astype(complex)
G1=np.zeros((4,4),complex); G2=np.zeros((4,4),complex)
for a,b in [(0,3),(1,2),(2,1),(3,0)]: G1[a,b]=1
for a,b,c in [(0,2,-1),(1,3,1),(2,0,-1),(3,1,1)]: G2[a,b]=c
G5=np.block([[np.zeros((2,2)),np.eye(2)],[np.eye(2),np.zeros((2,2))]]).astype(complex)
assert np.max(np.abs(G1@G1-I))<1e-15
assert np.max(np.abs(G2@G2-I))<1e-15
assert np.max(np.abs(G1@G2+G2@G1))<1e-15
assert np.max(np.abs(gamma0@G1+G1@gamma0))<1e-15
assert np.max(np.abs(gamma0@G2+G2@gamma0))<1e-15
assert np.max(np.abs(G5@G1-G1@G5))<1e-15
assert np.max(np.abs(G5@G2-G2@G5))<1e-15
Qminus=(I-G2)/2

def positive(H):
    es, v=np.linalg.eigh(H)
    return v[:,2:], (v[:,2:]@v[:,2:].conj().T)

def massive(k,lam,m,r0=.9,rtol=2e-11):
    x0=math.sqrt(r0)
    a=math.sqrt((1-r0)/r0)
    H=(k/a)*G1+(lam/r0)*G2+m*gamma0
    A,P=positive(H)
    def rhs(x,y):
        mat=2*k*x**3/(1-x*x)*G1+2*lam/math.sqrt(1-x*x)*G2+2*m*x*x/math.sqrt(1-x*x)*gamma0
        return (1j*mat@y.reshape(4,2)).ravel()
    sol=solve_ivp(rhs,(x0,0),A.ravel(),method='DOP853',rtol=rtol,atol=rtol*1e-2)
    if not sol.success: raise RuntimeError(sol.message)
    end=sol.y[:,-1].reshape(4,2)
    C=end.conj().T@Qminus@end
    lam_eig=np.linalg.eigvalsh((C+C.conj().T)/2)
    return {'k':k,'lambda':lam,'mass':m,'r0':r0,'eigenvalues':lam_eig.tolist(), 'average':float(lam_eig.mean()),'splitting':float(lam_eig[-1]-lam_eig[0]),'norm_error':float(np.max(np.abs(end.conj().T@end-np.eye(2)))), 'nfev':sol.nfev}

def universal(Lambda,rho0=25,rtol=2e-11):
    x0=math.sqrt(rho0)
    H=Lambda*(rho0*G1+rho0**(-.5)*G2)
    A,P=positive(H)
    def rhs(x,y): return (2j*Lambda*(x**3*G1+G2)@y.reshape(4,2)).ravel()
    sol=solve_ivp(rhs,(x0,0),A.ravel(),method='DOP853',rtol=rtol,atol=rtol*1e-2)
    end=sol.y[:,-1].reshape(4,2)
    C=end.conj().T@Qminus@end
    return {'Lambda':Lambda,'rho0':rho0,'eigenvalues':np.linalg.eigvalsh((C+C.conj().T)/2).tolist(), 'nfev':sol.nfev}

# Verify exponent bookkeeping EXACTLY using rational SymPy exponents.
Lsym,Ksym=sp.symbols('L K', positive=True)
Esym=(Lsym**2/Ksym)**sp.Rational(1,3)
rsym=(Lsym/Ksym)**sp.Rational(2,3)
rhosym=Esym**(-sp.Rational(2,3))
Lamsym=Lsym**sp.Rational(4,3)*Ksym**(-sp.Rational(1,3))
identity_exponents={
    'r2':(rsym*rhosym,Lsym**sp.Rational(2,9)*Ksym**(-sp.Rational(4,9))),
    'ratio_mass_theta':(1/(Ksym*sp.sqrt(rsym*rhosym)),Lsym**(-sp.Rational(1,9))*Ksym**(-sp.Rational(7,9))),
    'ratio_mass_proj':(rsym**sp.Rational(3,2)*rhosym**(-sp.Rational(1,2))/Lamsym,Lsym**(-sp.Rational(1,9))*Ksym**(-sp.Rational(7,9))),
    'mass_Duhamel':(rsym**sp.Rational(3,2)*rhosym**sp.Rational(3,2),Lsym**sp.Rational(1,3)*Ksym**(-sp.Rational(2,3))),
    'eps_match':(Lamsym*rsym*rhosym**3,Esym),
    'eps_tail':(rhosym**(-sp.Rational(3,2)),Esym),
}
for nm,(a,b) in identity_exponents.items():
    assert sp.simplify(sp.powdenest(a,force=True)/sp.powdenest(b,force=True)-1)==0,nm

# Verify endpoint geometry algebra and proper asymptotic estimates for a broad numeric grid.
checks=0; worst_ratio=0.; worst_geom=0.
for r0 in [.2,.5,.9]:
    e0=min(1/8, (r0/2)**(.75))
    for mass in [0,.5,3]:
      for k in [1,10,100,1000,10000]:
        for lam in [1,2,5,10,30]:
          E=(lam**2/k)**(1/3)
          if E>=e0: continue
          rstar=(lam/k)**(2/3); rho2=E**(-2/3); r2=rstar*rho2
          assert r2<r0/2 and rho2>=4-1e-10
          term_mass=1/(k*math.sqrt(r2))
          assert term_mass<=E*(1+1e-12)
          matchmass=(rstar**1.5*rho2**(-.5))/(lam**(4/3)*k**(-1/3))
          assert matchmass<=E*(1+1e-12)
          integral_mass=rstar**1.5*rho2**1.5
          assert integral_mass<=E*(1+1e-12)
          worst_ratio=max(worst_ratio,term_mass/E)
          checks+=1

# Differentiated angle: total variation bound evaluated by finite mesh quadrature.
from scipy.integrate import quad
case={'k':1e6,'lam':1.,'mass':1.,'r0':.9}
k_angle,lam_angle,m_angle,r0_angle=case.values()
k,lam,m,r0=k_angle,lam_angle,m_angle,r0_angle
E=(lam**2/k)**(1/3)
r2=(lam/k)**(2/3)*E**(-2/3)
def W(r): return math.sqrt(k*k*r/(1-r)+lam*lam/(r*r))
def beta(r): return math.atan2(m,W(r))
mesh=np.geomspace(r2,r0,10001)
variation=float(np.sum(np.abs(np.diff([beta(r) for r in mesh]))))
bound=2*m/(k*math.sqrt(r2))
assert variation<=bound*(1+1e-8)
# The previous false bound: peak is near r2, not r0.
old_bound=math.pi/2 *math.sin(beta(r0))
assert variation>old_bound

univ=[universal(v, rho0=35) for v in [.1,.5,1.,2.]]
max_chiral_splitting=max(max(z['eigenvalues'])-min(z['eigenvalues']) for z in univ)
assert max_chiral_splitting<1e-9

massive_results=[]
for m in [0,.5,2]:
  for k in [100,300,1000]:
    massive_results.append(massive(k,1,m))

output={
 'symbolic_exponent_identities':list(identity_exponents.keys()),
 'numeric_parameter_grid_passes':checks,
 'angle_counterexample':{'k':k_angle,'lambda':lam_angle,'m':m_angle,'r0':r0_angle,'r2':r2,'total_variation_mass_angle':variation,'corrected_upper_bound':bound,'old_outer_endpoint_bound':old_bound},
 'universal_two_chiral_eigenvalues':univ,
 'maximum_universal_chiral_splitting':max_chiral_splitting,
 'massive_schwarzschild': massive_results,
 'status':'PASS'
}
(OUT/'massive_operator_diagnostics.json').write_text(json.dumps(output,indent=2))
print(json.dumps({k:v for k,v in output.items() if k!='massive_schwarzschild'},indent=2))
for z in massive_results: print(z)
