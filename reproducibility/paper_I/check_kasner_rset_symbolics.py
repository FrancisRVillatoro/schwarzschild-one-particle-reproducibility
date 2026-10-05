#!/usr/bin/env python3
"""Independent symbolic checks for the terminal Kasner RSET."""
import sympy as sp

u,L,m = sp.symbols("u L m", positive=True, real=True)  # m = mu^2
f0=L**2+m

def fall(a,n):
    z=sp.Integer(1)
    for j in range(n):
        z*=a-j
    return sp.simplify(z)

F=[f0]
for n in range(1,5):
    F.append(sp.expand(L**2*(u**2*fall(sp.Rational(2,3),n)
                     +(1-u**2)*fall(sp.Rational(-4,3),n))))
f,f1,f2,f3,f4=F
O=sp.sqrt(f)
O1=f1/(2*O)
O2=(2*f*f2-f1**2)/(4*f**sp.Rational(3,2))
O3=(4*f**2*f3-6*f*f1*f2+3*f1**3)/(8*f**sp.Rational(5,2))
O4=(-4*(4*f1*f3+3*f2**2)*f**2+8*f**3*f4
    +36*f*f1**2*f2-15*f1**4)/(16*f**sp.Rational(7,2))

S,S1,S2=sp.Rational(1,4),sp.Rational(-1,2),sp.Rational(3,2)
w2=(4*O**2*S-2*O*O2+3*O1**2)/(8*O**3)
w4=-(16*O**4*S**2+16*O**4*S2-48*O**3*S*O2
     -80*O**3*O1*S1-8*O**3*O4+152*O**2*S*O1**2
     +80*O**2*O1*O3+52*O**2*O2**2
     -396*O*O1**2*O2+297*O1**4)/(128*O**7)

A=sp.simplify(1/(2*O)-w2/(2*O**2)+w2**2/(2*O**3)-w4/(2*O**2))
J=sp.factor(sp.integrate((3*u**2-1)*A,(u,0,1)))
print("J_mu =",J)
print("J_0  =",sp.factor(J.subs(m,0)))
assert sp.simplify(J.subs(m,0)+1/(63*L**3)) == 0

x=sp.symbols("x", positive=True)
I=sp.integrate(x**6*(224*x**4+1498*x**2+12089)
               /(1+x**2)**sp.Rational(13,2),(x,0,sp.oo))
res=sp.simplify(-I/(4*sp.pi**2*45360))
print("IR integral =",I)
print("Delta_IR   =",res)
assert sp.simplify(res + sp.Rational(523,498960)/sp.pi**2) == 0
