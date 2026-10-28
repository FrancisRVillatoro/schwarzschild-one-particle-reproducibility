#!/usr/bin/env python3
"""Exact algebra used in Work3. These finite checks supplement the proofs."""
import json
from pathlib import Path
import sympy as s
O=s.zeros(2); I=s.eye(2)
s1=s.Matrix([[0,1],[1,0]]);s2=s.Matrix([[0,-s.I],[s.I,0]]);s3=s.diag(1,-1)
g0=s.diag(1,1,-1,-1)
g=[g0]+[O.row_join(v).col_join((-v).row_join(O)) for v in [s1,s2,s3]]
G1=g0*g[1];G2=-g0*g[3];S1=s.diag(s1,s1);G5=O.row_join(I).col_join(I.row_join(O));I4=s.eye(4)
checks=[]
def equal(name,a,b):
 d=a-b
 ok=all(s.simplify(v)==0 for v in d) if isinstance(d,s.MatrixBase) else s.simplify(d)==0
 if not ok:raise AssertionError((name,d))
 checks.append(name)
x=s.symbols('x',real=True)
for l in range(1,9):
 p=s.jacobi(l-1,0,1,x);q=s.jacobi(l-1,1,0,x)
 equal(f'Jacobi pair 1 lambda={l}',p+(1+x)*s.diff(p,x),l*q)
 equal(f'Jacobi pair 2 lambda={l}',q-(1-x)*s.diff(q,x),l*p)
lam=s.symbols('lambda',integer=True,positive=True)
for n in range(1,7):
 equal(f'Hypergeometric centered product n={n}',s.rf(1-lam,n)*s.rf(1+lam,n),(-1)**n*s.prod(lam**2-j**2 for j in range(1,n+1)))
s0,d=s.symbols('s0 d',real=True);GG=s0*I4+d*S1
anti=lambda a,b:a*b+b*a
comm=lambda a,b:a*b-b*a
equal('G commutes G1',comm(GG,G1),s.zeros(4))
equal('G commutes g0',comm(GG,g0),s.zeros(4))
equal('Transverse Gram anticommutator',anti(GG,G2),2*s0*G2)
equal('Longitudinal Gamma identity',S1*G1,G5)
equal('Flux angular anticommutator',anti(GG*G1,G2),2*d*G5*G2)
equal('Flux mass anticommutator',anti(GG*G1,g0),s.zeros(4))
a,b,R,Rs=s.symbols('a b R Rs',positive=True)
k,kp,m,l=s.symbols('k kp mu lam',real=True)
H=lambda kk:kk/a*G1+l/R*G2+m*g0
K=(GG*H(kp)+H(k)*GG)/2
Kexpected=s0*H((k+kp)/2)+d*((k+kp)/(2*a)*G5+m*g0*S1)
equal('Full continuum energy kernel',K,Kexpected)
AA=[I4]+[-g0*gi for gi in g[1:]]
bs=s.symbols('b1:4',real=True)
C=[None]+[-bs[i-1]*g0*g[i]/2 for i in [1,2,3]]
for i,j in [(1,2),(1,3),(2,3)]:
 mat=s.I/4*(AA[i]*C[j]-C[j].conjugate().T*AA[i]+AA[j]*C[i]-C[i].conjugate().T*AA[j])
 equal(f'Shear connection ({i},{j})',mat,(bs[i-1]-bs[j-1])*s.I*g[i]*g[j]/4)
a2=(Rs-R)/R
coef=(-Rs/(2*R**2)-a2/R)/(4*a2*R**2)
equal('Physical shear coefficient',coef,-(3*Rs-2*R)/(8*R**3*(Rs-R)))
equal('Shear endpoint coefficient',s.limit(R**3*coef,R,0),-s.Rational(3,8))
Ip=-R-Rs*s.log(1-R/Rs)
equal('Longitudinal comparison primitive',s.diff(Ip,R),1/a2)
equal('Comparison leading coefficient',s.limit(Ip/R**2,R,0),1/(2*Rs))
# length exponents: dκ=-1, f=+1/2, B=-1/2, q=aR²=+2.
equal('Field length exponent',-s.Integer(1)-s.Rational(1,2),-s.Rational(3,2))
equal('Energy length exponent',-s.Rational(1,2)-s.Rational(5,2)-1,-4)
equal('Shear length exponent',-3-1,-4)
summary={'exact_checks_passed':len(checks),'checks':checks}
out=Path(__file__).with_name('work3_symbolic_results.json');out.write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
