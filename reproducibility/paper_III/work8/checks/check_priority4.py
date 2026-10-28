#!/usr/bin/env python3
"""Independent algebraic audit: CAR global/local purity and rank-one excitations.
This script checks exact identities, not Hadamard microlocal estimates.
"""
import json
import numpy as np
import sympy as s
from pathlib import Path

base=Path(__file__).resolve().parents[1]
results={}

def check(name,test,value=None):
    ok=bool(test)
    results[name]={'pass':ok}
    if value is not None:results[name]['value']=str(value)
    if not ok:raise AssertionError(f'FAILED: {name}: {value}')

sq2=s.sqrt(2);sq3=s.sqrt(3)
u1=s.Matrix([1,0,1,0])/sq2
u2=s.Matrix([0,sq3/s.Integer(2),0,s.Rational(1,2)])
v1=s.Matrix([1,0,-1,0])/sq2
J=s.Matrix([[1,0],[0,1],[0,0],[0,0]])
P=u1*u1.T+u2*u2.T
E=J*J.T
C=J.T*P*J
Pp=P-u1*u1.T
Pm=P+v1*v1.T
Cp=J.T*Pp*J
Cm=J.T*Pm*J
Z4=s.zeros(4)
Z2=s.zeros(2)
I4=s.eye(4)
I2=s.eye(2)
check('global_pure_projection_exact',s.simplify(P*P-P)==Z4)
check('isometric_embedding',J.T*J==I2)
check('global_perturbations_are_projections', Pp*Pp==Pp and Pm*Pm==Pm)
check('orthogonality_particle_antiparticle',u1.dot(u2)==0 and P*v1==s.zeros(4,1))
check('restricted_covariance_exact',C==s.diag(s.Rational(1,2),s.Rational(3,4)))
check('purity_defect_factorization',s.simplify(C-C*C-J.T*P*(I4-E)*P*J)==Z2)
check('purity_defect_nonzero',C-C*C==s.diag(s.Rational(1,4),s.Rational(3,16)))
check('particle_covariance_compression',Cp==C-(J.T*u1)*(J.T*u1).T)
check('antiparticle_covariance_compression',Cm==C+(J.T*v1)*(J.T*v1).T)
check('particle_excited_local_mixed',Cp==s.diag(0,s.Rational(3,4)))
check('antiparticle_excited_local_mixed',Cm==s.diag(1,s.Rational(3,4)))
check('local_fermion_covariance_bounds',all(0<=d<=1 for M in (C,Cp,Cm) for d in M.diagonal()))
check('restriction_fraction', (J.T*u1).dot(J.T*u1)==s.Rational(1,2))

# Direct CAR Fock calculation. Two particle modes, two antiparticle modes;
# Dirac field components = (a1,a2,b1^dagger,b2^dagger).
Id=np.eye(2,dtype=int)
Z=np.diag([1,-1])
ann=np.array([[0,1],[0,0]],dtype=int)
a=[]
for k in range(4):
    o=np.array([[1]],dtype=int)
    for j in range(4):o=np.kron(o,Z if j<k else ann if j==k else Id)
    a.append(o)
create=[t.T for t in a]
fullId=np.eye(16,dtype=int)
for i in range(4):
    for j in range(4):
        lhs=a[i]@create[j]+create[j]@a[i]
        check(f'CAR_{i}_{j}',np.array_equal(lhs,fullId if i==j else np.zeros((16,16),dtype=int)))

vac=np.zeros((16,1),dtype=int);vac[0,0]=1
psi=[a[0],a[1],create[2],create[3]]

def cov(ket):
    return np.array([[(ket.T@psi[i]@psi[j].T@ket)[0,0] for j in range(4)] for i in range(4)],dtype=int)

k_p=create[0]@vac
k_a=create[2]@vac
check('vacuum_norm',int((vac.T@vac)[0,0])==1)
check('particle_norm',int((k_p.T@k_p)[0,0])==1)
check('antiparticle_norm',int((k_a.T@k_a)[0,0])==1)
P0=np.diag([1,1,0,0])
Ppart=np.diag([0,1,0,0])
Panti=np.diag([1,1,1,0])
check('Fock_vacuum_CAR_particle_covariance',np.array_equal(cov(vac),P0))
check('Fock_excited_particle_subtracts_rank_one',np.array_equal(cov(k_p),Ppart))
check('Fock_excited_antiparticle_adds_rank_one',np.array_equal(cov(k_a),Panti))

# Necessary compactness: sqrt(C)-sqrt(D) Hilbert-Schmidt implies C-D HS.
# Check identity on noncommuting 2x2 contractions.
from scipy.linalg import sqrtm
C1=np.array([[0.8,0.1],[0.1,0.3]])
C2=np.array([[0.4,-0.17],[-0.17,0.6]])
rC1=sqrtm(C1);rC2=sqrtm(C2)
d=rC1-rC2
prod=d@rC1+rC2@d
check('sqrt_HS_ideal_identity_noncommuting',np.linalg.norm(C1-C2-prod)<1.e-12, np.linalg.norm(C1-C2-prod))

results['total_checks']=len(results)
results['all_pass']=all(v.get('pass',True) for k,v in results.items() if isinstance(v,dict))
(base/'results'/'priority4_exact_CAR_checks.json').write_text(json.dumps(results,indent=2,ensure_ascii=False))
print('EXACT CAR CHECKS:',results['total_checks'],'ALL PASS:',results['all_pass'])
print('Global P rank',P.rank(),'; restricted eigenvalues',list(C.diagonal()))
print('Local purity defect',list((C-C*C).diagonal()))
print('Local particle excited eigenvalues',list(Cp.diagonal()))
print('Local antiparticle excited eigenvalues',list(Cm.diagonal()))
