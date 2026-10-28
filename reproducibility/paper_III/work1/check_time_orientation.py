#!/usr/bin/env python3
"""Independent Work-1 algebra and evolution checks. Does not test the missing
spectral-symbol estimates or the full physical Unruh-packet construction.
Requires numpy, scipy, sympy. Run: python check_time_orientation.py
Outputs JSON and human-readable diagnostics beside this script.
"""
from __future__ import annotations
import json
import platform
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import expm
import sympy as sp

HERE = Path(__file__).resolve().parent
report: dict = {"environment": {"python": platform.python_version(), "numpy": np.__version__,
                               "scipy": scipy.__version__, "sympy": sp.__version__},
                "symbolic": {}, "mode_tests": [], "universal_tests": []}
I = sp.I
sig = [sp.Matrix([[0,1],[1,0]]), sp.Matrix([[0,-I],[I,0]]), sp.diag(1,-1)]
Z = sp.zeros(2)
beta = sp.diag(1,1,-1,-1)
gamma = [beta] + [sp.BlockMatrix([[Z,s],[-s,Z]]).as_explicit() for s in sig]
g1 = beta*gamma[1]
g2 = -beta*gamma[3]
g5 = sp.BlockMatrix([[Z,sp.eye(2)],[sp.eye(2),Z]]).as_explicit()
s1 = sp.diag(sig[0],sig[0]); s3 = sp.diag(sig[2],sig[2])
comm = lambda A,B: A*B-B*A
anti = lambda A,B: A*B+B*A

def zero(name: str, expression):
    if isinstance(expression, sp.MatrixBase):
        ok = all(sp.simplify(x) == 0 for x in expression)
    else:
        ok = sp.simplify(expression) == 0
    report["symbolic"][name] = bool(ok)
    if not ok:
        raise AssertionError(f"Symbolic check failed: {name}: {expression}")

p,q,m,pp,s,d,E,v = sp.symbols('p q m pp s d E v', real=True)
H = p*g1+q*g2+m*beta
Hp = pp*g1+q*g2+m*beta
zero('H_squared',H*H-(p*p+q*q+m*m)*sp.eye(4))
zero('lambda_flip_mass_fixed',s1*H*s1-(p*g1-q*g2+m*beta))
zero('k_flip_mass_fixed',s3*H*s3-(-p*g1+q*g2+m*beta))
zero('Gamma1_also_flips_mass',g1*H*g1-(p*g1-q*g2-m*beta))
zero('chirality_commutes_massless_H',comm(g5,p*g1+q*g2))
zero('chirality_anticommutes_beta',anti(g5,beta))
P0 = (sp.eye(2)+sig[2])/2
P0p = v*sig[0]/2
Xminus = -I/(2*E)*comm(P0,P0p)
zero('physical_first_superadiabatic_correction',Xminus-v*sig[1]/(4*E))
zero('physical_first_order_transport',-I*P0p-comm(E*sig[2],Xminus))
for S in sig[:2]:
    zero('J_commutes_sigma_'+str(sig.index(S)+1),sig[0]*S.conjugate()*sig[0]-S)

# Spin connection derived from tetrad and Christoffel symbols, not inserted.
tau,t,theta,phi = sp.symbols('sigma t theta phi',real=True)
a=sp.Function('a')(tau); r=sp.Function('r')(tau)
coords=[tau,t,theta,phi]
g=sp.diag(1,-a*a,-r*r,-r*r*sp.sin(theta)**2)
ginv=sp.simplify(g.inv())
coframe=sp.diag(-1,a,r,r*sp.sin(theta)); invframe=coframe.inv()
eta=sp.diag(1,-1,-1,-1)
Christ={}
for A in range(4):
 for mu in range(4):
  for nu in range(4):
   Christ[A,mu,nu]=sp.simplify(sum(ginv[A,k]*(sp.diff(g[k,nu],coords[mu])+sp.diff(g[k,mu],coords[nu])-sp.diff(g[mu,nu],coords[k]))/2 for k in range(4)))
spin=[]
for mu in range(4):
 Om=sp.zeros(4)
 for A in range(4):
  for B in range(4):
   Om[A,B]=sp.simplify(sum(coframe[A,nu]*(sp.diff(invframe[nu,B],coords[mu])+sum(Christ[nu,mu,lam]*invframe[lam,B] for lam in range(4))) for nu in range(4)))
 lower=eta*Om
 G=sp.zeros(4)
 for A in range(4):
  for B in range(4):
   G+=lower[A,B]*gamma[A]*gamma[B]/4
 spin.append(sp.simplify(G))
zero('spin_connection_sigma',spin[0])
zero('spin_connection_t',spin[1]+sp.diff(a,tau)/2*gamma[0]*gamma[1])
zero('spin_connection_theta',spin[2]+sp.diff(r,tau)/2*gamma[0]*gamma[2])
zero('spin_connection_phi',spin[3]+sp.diff(r,tau)*sp.sin(theta)/2*gamma[0]*gamma[3]-sp.cos(theta)/2*gamma[2]*gamma[3])
contract=sp.zeros(4)
for mu in range(4):
 for A in range(4):
  contract+=gamma[A]*invframe[mu,A]*spin[mu]
expected=-gamma[0]*(sp.diff(a,tau)/(2*a)+sp.diff(r,tau)/r)+gamma[2]*sp.cot(theta)/(2*r)
zero('contracted_spin_connection',contract-expected)

# Exact local bilinear identities, retaining the angular Gram matrix.
G=s*sp.eye(4)+d*s1
pbar=(p+pp)/2
Hbar=pbar*g1+q*g2+m*beta
T00=(G*Hp+H*G)/2
T11=pbar*G*g1
T01=-(p+pp)*G/4-(G*g1*Hp+H*G*g1)/4
zero('T00_local_bilinear',T00-(s*Hbar+d*(pbar*g5+m*beta*s1)))
zero('T11_local_bilinear',T11-pbar*(s*g1+d*g5))
zero('T01_local_bilinear',T01-(-pbar*G-q*d*g5*g2/2))
zero('T01_transverse_anticommutator',anti(G*g1,g2)-2*d*g5*g2)
zero('T01_mass_anticommutator',anti(G*g1,beta))
report['symbolic']['T01_transverse_term_not_zero'] = bool(anti(G*g1,g2) != sp.zeros(4))

# Static calibration: positive physical energy and the remaining-time phase.
HH=np.array((sp.Rational(7,10)*g1+sp.Rational(9,10)*g2+sp.Rational(1,5)*beta),complex)
eig,mat=np.linalg.eigh(HH); omega=float(eig[-1]); pplus=(np.eye(4)+HH/omega)/2
pos=mat[:,-1]; neg=mat[:,0]
def static_energy(z):
    d_sigma=1j*HH@z
    return float((-1j/2*(np.vdot(z,d_sigma)-np.vdot(d_sigma,z))).real)
CC=expm(1j*HH*.47)@pplus
report['static_calibration']={'positive_particle_energy':static_energy(pos),'expected_omega':omega,
    'negative_mode_classical_energy':static_energy(neg),'antiparticle_energy':-static_energy(neg),
    'two_point_equation_error':float(np.linalg.norm(-1j*(1j*HH@CC)-HH@CC)),
    'positive_projection_error':float(np.linalg.norm(pplus@pplus-pplus))}
assert abs(static_energy(pos)-omega)<1e-12
assert abs(static_energy(neg)+omega)<1e-12

# Finite CAR matrix test of rank-one signs.
paulis=[np.array(x,dtype=complex) for x in sig]
ann=np.array([[0.,1.],[0.,0.]])
aa=np.kron(ann,np.eye(2)); bb=np.kron(np.diag([1.,-1.]),ann)
vac=np.array([1.,0.,0.,0.],complex)
part=aa.conj().T@vac; apart=bb.conj().T@vac
rng=np.random.default_rng(49017)
uvec=rng.normal(size=4)+1j*rng.normal(size=4)
vvec=rng.normal(size=4)+1j*rng.normal(size=4)
ops=[uvec[A]*aa+vvec[A]*bb.conj().T for A in range(4)]
def car_cov(state):
 return np.array([[np.vdot(state,ops[A]@ops[B].conj().T@state) for B in range(4)] for A in range(4)])
report['CAR_particle_covariance_error']=float(np.max(np.abs(car_cov(part)-car_cov(vac)+np.outer(uvec,uvec.conj()))))
report['CAR_antiparticle_covariance_error']=float(np.max(np.abs(car_cov(apart)-car_cov(vac)-np.outer(vvec,vvec.conj()))))
assert report['CAR_particle_covariance_error']<1e-12
assert report['CAR_antiparticle_covariance_error']<1e-12

A1=np.array(g1,dtype=complex); A2=np.array(g2,dtype=complex)
B0=np.array(beta,dtype=complex); G5=np.array(g5,dtype=complex)
S1=np.array(s1,dtype=complex); Gphys=[np.array(x,dtype=complex) for x in gamma]
E4=np.eye(4,dtype=complex)
def h(r,k,lam,m):
 return k*np.sqrt(r/(1-r))*A1+lam/r*A2+m*B0

def evolve(k,lam,m,sgn):
 """i*sgn*dchi/dsigma=H chi. Work in z=sqrt(r), regular at z=0."""
 z0=np.sqrt(.9)
 def rhs(z,state):
  U=state.reshape(4,4)
  Az=2*k*z**3/(1-z*z)*A1+2*lam/np.sqrt(1-z*z)*A2+2*m*z*z/np.sqrt(1-z*z)*B0
  return (-1j*sgn*Az@U).ravel()
 sol=solve_ivp(rhs,(z0,0.),E4.ravel(),method='DOP853',rtol=2e-11,atol=2e-13)
 if not sol.success: raise RuntimeError(sol.message)
 return sol.y[:,-1].reshape(4,4)

for k,lam,mass in [(5,2,0),(30,3,0),(5,2,1),(100,1,0),(-5,2,1),(5,-2,1)]:
 U=evolve(k,lam,mass,-1); Ua=evolve(k,lam,mass,+1)
 vals,vecs=np.linalg.eigh(h(.9,k,lam,mass))
 C=vecs[:,vals>0]
 Qm=(E4-np.sign(lam)*A2)/2
 probs=[float(np.linalg.norm(Qm@U@C[:,j])**2) for j in range(2)]
 probsaux=[float(np.linalg.norm(Qm@Ua@C[:,j].conj())**2) for j in range(2)]
 P=C@C.conj().T; Pterminal=U.conj().T@(E4-Qm)@U
 errtype=np.linalg.norm(U@(E4-Pterminal)@P-Qm@U@P)
 row={'k':k,'lambda':lam,'mass':mass,'physical_occupations':probs,
      'auxiliary_occupations':probsaux,'max_occupation_change':max(abs(x-y) for x,y in zip(probs,probsaux)),
      'conjugate_propagator_error':float(np.linalg.norm(U-Ua.conj())),
      'unitarity_error':float(np.linalg.norm(U.conj().T@U-E4)),
      'typed_Bogoliubov_identity_error':float(errtype)}
 assert row['max_occupation_change']<1e-10
 assert row['conjugate_propagator_error']<1e-10
 assert row['unitarity_error']<2e-8
 assert row['typed_Bogoliubov_identity_error']<2e-8
 report['mode_tests'].append(row)
 print('MODE',row,flush=True)

# Universal finite-endpoint regression: not a precision recalculation of F(infinity).
for Lam in [.1,.5,1.]:
 eta0=2*np.sqrt(25.)
 h0=Lam*(paulis[1]+eta0**3/8*paulis[0])
 vals,vecs=np.linalg.eigh(h0); aux0=vecs[:,1]; phys0=paulis[0]@aux0.conj()
 def rhs(eta,z):
  H=Lam*(paulis[1]+eta**3/8*paulis[0])
  return np.concatenate([1j*H@z[:2],-1j*H@z[2:]])
 sol=solve_ivp(rhs,(eta0,0),np.r_[phys0,aux0],method='DOP853',rtol=2e-11,atol=2e-13)
 if not sol.success:raise RuntimeError(sol.message)
 ph,aux=sol.y[:2,-1],sol.y[2:,-1]
 Qm=(np.eye(2)-paulis[1])/2
 row={'Lambda':Lam,'rho_start':25.,'physical_probability':float(np.linalg.norm(Qm@ph)**2),
      'auxiliary_probability':float(np.linalg.norm(Qm@aux)**2),
      'J_mapping_error':float(np.linalg.norm(ph-paulis[0]@aux.conj()))}
 assert row['J_mapping_error']<1e-10
 report['universal_tests'].append(row)
 print('UNIVERSAL',row,flush=True)

# Random complex packet, direct Dirac tensor vs reduced formulas.
maxerr=0.
for trial in range(12):
 theta0=rng.normal(size=2)+1j*rng.normal(size=2)
 theta0/=np.linalg.norm(theta0)
 T=np.column_stack([theta0,paulis[0]@theta0])
 M=np.kron(np.eye(2),T)
 st=float(np.vdot(theta0,theta0).real); dt=float(np.vdot(theta0,paulis[0]@theta0).real)
 r0=.19;a0=np.sqrt(1/r0-1); mass=.3;lam=2.;adot=-.7
 ks=np.array([2.,3.5]);co=rng.normal(size=2)+1j*rng.normal(size=2)
 chi=rng.normal(size=(2,4))+1j*rng.normal(size=(2,4))
 ph=np.exp(1j*ks*.27); ds=co*ph
 W=sum(ds[j]*(M@chi[j]) for j in range(2))
 der0=sum(ds[j]*(-1j*M@h(r0,ks[j],lam,mass)@chi[j]) for j in range(2))+.123*W
 spin1=-adot/(2*a0)*Gphys[0]@Gphys[1]
 der1=sum(ds[j]*1j*ks[j]/a0*(M@chi[j]) for j in range(2))+spin1@W
 ders=[der0,der1]; lower=[Gphys[0],-Gphys[1]]
 direct=[]
 for A,B in [(0,0),(1,1),(0,1)]:
  value=1j/4*(np.vdot(W,B0@lower[A]@ders[B])+np.vdot(W,B0@lower[B]@ders[A])
              -np.vdot(ders[A],B0@lower[B]@W)-np.vdot(ders[B],B0@lower[A]@W))
  direct.append(value.real)
 reduced=np.zeros(3,dtype=complex)
 for j in range(2):
  for jj in range(2):
   kb=(ks[j]+ks[jj])/2
   HM=h(r0,kb,lam,mass)
   mats=[st*HM+dt*(kb/a0*G5+mass*B0@S1),kb/a0*(st*A1+dt*G5),
         -kb/a0*(st*E4+dt*S1)-lam*dt/(2*r0)*G5@A2]
   for n in range(3):
    reduced[n]+=ds[j].conjugate()*ds[jj]*np.vdot(chi[j],mats[n]@chi[jj])
 err=float(np.max(np.abs(np.array(direct)-reduced.real)))
 maxerr=max(maxerr,err)
 assert err<1e-10
report['random_local_stress_absolute_error']=maxerr

# Explicit non-cancellation counterexample for T01.
chi=np.array([0,1,0,1],complex)/np.sqrt(2)
sval=1.;dval=np.sqrt(3)/2;rval=.1;aval=3.;lamval=1.
report['T01_counterexample']={'k':0,'lambda':1,'mass':0,'r':rval,'s_Theta':sval,'d_Theta':dval,
 'rho':float((lamval/(aval*rval**3)*np.vdot(chi,A2@chi)).real),
 'T01':float((-lamval*dval/(2*aval*rval**3)*np.vdot(chi,G5@A2@chi)).real),'T11':0.}
report['status']='all_implemented_checks_passed'
(HERE/'Work1_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print('SYMBOLIC CHECKS:',len(report['symbolic']), 'PASS',flush=True)
print('CAR rank-one signs: PASS',flush=True)
print('Random stress max absolute error:',maxerr,flush=True)
print('T01 counterexample:',report['T01_counterexample'],flush=True)
print('ALL IMPLEMENTED WORK-1 CHECKS PASSED. This is not a proof of the spectral-symbol estimates.',flush=True)
