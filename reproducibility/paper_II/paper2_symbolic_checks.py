"""Symbolic and high-precision checks for paper II (units 2M = 1).
 1. Pump Q = (sqrt q)''/sqrt q, friction H = q'/q, its zero at r = 3/4 (maximal-volume slice) and Q(3/4).
 2. Velocity-count Bogoliubov map (Proposition 3): |alpha_v|^2 - |beta_v|^2 = 1.
 3. Normalization of the comparison solutions: W[u,u*] = (2/pi)(P conj(Q) - conj(P) Q), hence 2 Im(conj(P) Q) = -pi/2.
 4. Weber-Schafheitlin integral, the constants D and A_s (Proposition 6), and the Abel-regularized integral of z H0^2.
 5. Fixed-l matching algebra (Theorem 7): alpha - beta = 2, alpha + beta, |alpha|^2 - |beta|^2 = 1, |beta|^2.
 6. Curvature-coupling term of the excitation stress (Remark 13) on the Kasner limit."""
import sympy as sp, mpmath as mp, numpy as np
# 1
r=sp.symbols('r',positive=True); f=1/r-1; q=r**2*sp.sqrt(f); W=sp.sqrt(q)
d=lambda F: sp.sqrt(f)*sp.diff(F,r)
Q=sp.simplify(d(d(W))/W); H=sp.simplify(d(q)/q)
print('Q =',sp.factor(Q),'  Q - formula (2.4):',sp.simplify(Q+(1/r**3)*(sp.Rational(1,2)+1/(16*(1-r)))))
print('H =',sp.factor(H),'  zeros:',sp.solve(sp.numer(sp.together(H)),r),'  Q(3/4) =',Q.subs(r,sp.Rational(3,4)))
print('Q = H_dot/2 + H^2/4 check:',sp.simplify(Q-(d(H)/2+H**2/4)))
# 2
om,Om,Hs=sp.symbols('omega Omega H',positive=True); s_=sp.sqrt(Om/om)
al=(s_+1/s_)/2+sp.I*Hs/(4*sp.sqrt(om*Om)); be=(s_-1/s_)/2+sp.I*Hs/(4*sp.sqrt(om*Om))
print('velocity map: |alpha_v|^2-|beta_v|^2 =',sp.simplify(sp.expand(al*sp.conjugate(al)-be*sp.conjugate(be))))
# 3
mp.mp.dps=20
P,Qc,z=mp.mpc(1,2),mp.mpc(0.5,-1),mp.mpf('1.3')
u=lambda zz: mp.sqrt(zz)*(P*mp.besselj(0,zz)+Qc*mp.bessely(0,zz))
Wr=u(z)*mp.conj(mp.diff(u,z))-mp.diff(u,z)*mp.conj(u(z))
print('W[u,u*] =',mp.nstr(Wr,12),'  (2/pi)(P Qbar - Pbar Q) =',mp.nstr(2/mp.pi*(P*mp.conj(Qc)-mp.conj(P)*Qc),12),'  -> W=i requires 2 Im(Pbar Q) =',mp.nstr(-mp.pi/2,8))
# 4
fWS=lambda x: x**-0.5*(mp.besselj(0,x)-1j*mp.bessely(0,x))**2
I_WS=mp.quad(fWS,[0,1])+mp.quadosc(fWS,[1,mp.inf],period=mp.pi)
ref=-mp.sqrt(2*mp.pi)*mp.gamma(0.25)/mp.gamma(0.75)**3*(1-1j)
D=mp.pi**1.5*mp.gamma(0.25)/(24*mp.gamma(0.75)**3)
print('WS integral =',mp.nstr(I_WS,12),' closed form',mp.nstr(ref,12),'  D =',mp.nstr(D,12),'  A_s = 2D^2 =',mp.nstr(2*D**2,12))
from scipy.integrate import quad as squad
from scipy.special import hankel2
def fA(eps):
    re=squad(lambda x: (x*hankel2(0,x)**2*np.exp(-eps*x)).real,0,60/eps,limit=4000)[0]
    im=squad(lambda x: (x*hankel2(0,x)**2*np.exp(-eps*x)).imag,0,60/eps,limit=4000)[0]
    return re+1j*im
import numpy as np
vals=[fA(e) for e in (0.08,0.04,0.02,0.01)]
print('Abel-regularized int z H0(2)^2 e^{-eps z} dz, eps=0.08,0.04,0.02,0.01 ->',[f'{v.real:.6f}{v.imag:+.6f}i' for v in vals],'  expected 2/pi^2 =',f'{2/np.pi**2:.6f}')
# 5
K,L=sp.symbols('K L',positive=True); g=sp.EulerGamma
apb=sp.Rational(1,2)+sp.I/sp.pi*(4*sp.log(L)-sp.log(K/4)+3*g); amb=2
alH=(apb+amb)/2; beH=(apb-amb)/2; Lam=K/(4*sp.exp(3*g)*L**4)
print('alpha_H =',sp.simplify(alH-(sp.Rational(5,4)-sp.I/(2*sp.pi)*sp.log(Lam))),' (should be 0);  beta_H check:',sp.simplify(beH-(-sp.Rational(3,4)-sp.I/(2*sp.pi)*sp.log(Lam))))
print('|alpha|^2-|beta|^2 =',sp.simplify(sp.expand(alH*sp.conjugate(alH)-beH*sp.conjugate(beH))),'  |beta|^2 - [9/16 + ln^2 Lambda/(4pi^2)] =',sp.simplify(sp.expand(beH*sp.conjugate(beH))-(sp.Rational(9,16)+sp.log(Lam)**2/(4*sp.pi**2))))
# 6
y,X,Z1,Z2,a,b=sp.symbols('y X Z1 Z2 a b',positive=True); x=[y,X,Z1,Z2]
gK=sp.diag(-1,y**sp.Rational(-2,3),y**sp.Rational(4,3),y**sp.Rational(4,3)); gi=gK.inv()
G=[[[sum(gi[i,dd]*(sp.diff(gK[dd,j],x[c])+sp.diff(gK[dd,c],x[j])-sp.diff(gK[j,c],x[dd])) for dd in range(4))/2 for c in range(4)] for j in range(4)] for i in range(4)]
F=(a+b*sp.log(y))**2
Hm=sp.simplify(gi*sp.Matrix(4,4,lambda j,c: sp.diff(F,x[j],x[c])-sum(G[i][j][c]*sp.diff(F,x[i]) for i in range(4)))); box=sp.simplify(Hm.trace())
sh=[sp.expand(sp.simplify((box-Hm[i,i])*y**2)) for i in range(4)]
print('excitation coupling term, coefficient of b^2 ln y / y^2 in (delta Box - nabla nabla)F, components (y,X,Z,Z):',[sp.Poly(e,sp.log(y)).coeff_monomial(sp.log(y))/b**2 for e in sh])
