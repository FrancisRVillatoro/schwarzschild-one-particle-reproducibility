#!/usr/bin/env python3
"""Exact algebraic checks for Work 2 (not a numerical proof of symbol bounds).

Requires SymPy.  Run: python check_symbol_algebra.py
Checks the physical sign, the first two formal coefficients, their homogeneity,
the transverse-axis derivative, and the exponent bookkeeping of the proof.
No Schwarzschild campaign or finite-difference spectral differentiation is run.
"""
from __future__ import annotations
from fractions import Fraction
import sys
import sympy as sp

count = 0

def simp(x):
    return sp.factor(sp.simplify(x))

def check(name, expressions):
    global count
    if isinstance(expressions, sp.MatrixBase):
        expressions = list(expressions)
    elif not isinstance(expressions, (tuple, list)):
        expressions = [expressions]
    remainders = [simp(x) for x in expressions]
    bad = [x for x in remainders if x != 0]
    if bad:
        raise AssertionError(f'{name}: nonzero residuals {bad}')
    count += 1
    print(f'PASS {count:02d}: {name}')

z = sp.symbols('z', positive=True)   # y = z^3
P, Q1, Q2 = sp.symbols('P Q1 Q2', real=True)
R = sp.symbols('R', positive=True)
I = sp.I

def dy(expr):
    return expr.diff(z) / (3*z**2)

def euler(expr):
    return P*expr.diff(P) + Q1*expr.diff(Q1) + Q2*expr.diff(Q2)

d = sp.Matrix([P*z, Q1/z**2, Q2/z**2])
W2 = P**2*z**6 + Q1**2 + Q2**2
omega = sp.sqrt(W2)/z**2
n0 = (d/omega).applyfunc(simp)
# Physical covariance equation: n' = -2 d x n.
n1 = (d.cross(dy(n0))/(2*omega**2)).applyfunc(simp)
n2 = (d.cross(dy(n1))/(2*omega**2) - n0*(n1.dot(n1))/2).applyfunc(simp)

check('principal unit vector', n0.dot(n0)-1)
check('physical order-zero transport: n0_prime + 2 d cross n1', dy(n0)+2*d.cross(n1))
check('physical order-one transport: n1_prime + 2 d cross n2', dy(n1)+2*d.cross(n2))
check('projection constraint at order -1', n0.dot(n1))
check('projection constraint at order -2', 2*n0.dot(n2)+n1.dot(n1))
check('Euler homogeneity of n0 (degree 0)', euler(n0))
check('Euler homogeneity of n1 (degree -1)', euler(n1)+n1)
check('Euler homogeneity of n2 (degree -2)', euler(n2)+2*n2)

s1 = sp.Matrix([[0,1],[1,0]])
s2 = sp.Matrix([[0,-I],[I,0]])
s3 = sp.diag(1,-1)
zero = sp.zeros(2)
G = [sp.Matrix.vstack(sp.Matrix.hstack(zero,s),sp.Matrix.hstack(s,zero))
     for s in (s1,s2,s3)]
for i in range(3):
    check(f'Clifford involution {i+1}', G[i]*G[i]-sp.eye(4))
for i in range(3):
    for j in range(i+1,3):
        check(f'Clifford anticommutator {i+1},{j+1}', G[i]*G[j]+G[j]*G[i])
H = sum((d[i]*G[i] for i in range(3)), sp.zeros(4))
Pi0 = (sp.eye(4)+H/omega)/2
comm = lambda a,b:a*b-b*a
pi1 = (-I*comm(Pi0,dy(Pi0))/(2*omega)).applyfunc(simp)
pi1_explicit = I*P*z**2*G[0]*(Q1*G[1]+Q2*G[2])/(4*W2**sp.Rational(3,2))
check('explicit physical first homogeneous coefficient', pi1-pi1_explicit)
check('first formal transport in four-component representation', comm(H,pi1)+I*dy(Pi0))
check('first coefficient is tangent to the projection manifold', Pi0*pi1+pi1*Pi0-pi1)
check('transverse-axis derivative is y Gamma_1/(2R)',
      Pi0.diff(P).subs({P:0,Q1:R,Q2:0})-z**3*G[0]/(2*R))
check('first coefficient vanishes on the transverse axis', pi1.subs(P,0))
check('first coefficient vanishes on the longitudinal axis', pi1.subs({Q1:0,Q2:0}))

# Exact finite-endpoint identity in the physical sign: no numerical integration.
# d_s[U(y,s) A(s) U(s,y)] = U(A' - i[H,A])U = i U R_A U.
hx,hy,hz,a0,ax,ay,az,ap0,apx,apy,apz = sp.symbols(
    'hx hy hz a0 ax ay az ap0 apx apy apz', real=True)
HH = hx*s1+hy*s2+hz*s3
AA = a0*sp.eye(2)+ax*s1+ay*s2+az*s3
Aprime = ap0*sp.eye(2)+apx*s1+apy*s2+apz*s3
residual = -I*Aprime-comm(HH,AA)
check('sign of the conjugated residual identity', Aprime-I*comm(HH,AA)-I*residual)

n_cases = 0
n_terms = 0
for N in range(7):
    for p in range(5):
        for q in range(9):
            M = max(4*q+1, N+1+q+p)
            assert M > 4*q
            assert p-M <= -N-1-q
            worst = -1-Fraction(M-4*q,3)
            assert worst < -1
            for g in range(q+1):
                term = -1-Fraction(M,3)+g+Fraction(4*(q-g),3)
                assert term <= worst
                assert p-M-g <= p-M
                n_terms += 1
            n_cases += 1
print(f'PASS bookkeeping: {n_cases} (N,p,|alpha|) cases and {n_terms} derivative allocations')
print(f'Exact algebraic identities checked: {count}')
print('ALL CHECKS PASSED')
print('Scope: algebra and exponent bookkeeping only; the all-order proof is in the manuscript.')
print('Python:', sys.version.split()[0], 'SymPy:', sp.__version__)
