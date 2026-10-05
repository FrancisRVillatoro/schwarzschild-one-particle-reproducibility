"""
Equal-time (spatial) point splitting on the vacuum Kasner geometry
    g = -dy^2 + y^(-2/3) dX^2 + y^(4/3) (dZ1^2 + dZ2^2),  base point y = 1.

Computes, for the massless minimally coupled scalar,
  (i)  the finite part of the Hadamard parametrix H minus the Fourier transform of the
       massless second-order adiabatic subtraction, for <Phi^2>;
  (ii) the same for the anisotropic combination <T^X_X - T^Z_Z> with fourth-order
       subtraction.
The Hadamard coefficients are obtained here from their defining equations
(Hamilton-Jacobi equation for sigma, transport equations for U and V), not from tables.
Conventions: H = (1/8 pi^2) [U/sigma + V ln sigma],  V = V0 + V1 sigma + ...
"""
import sympy as sp

s, a, b, c, n = sp.symbols('s a b c n', real=True)
X = [s, a, b, c]
DMAX = 6

def deg_part(e, k):
    p = sp.Poly(sp.expand(e), *X)
    return sp.Add(*[co*s**m[0]*a**m[1]*b**m[2]*c**m[3] for m, co in p.terms() if sum(m) == k])

def trunc(e, k):
    p = sp.Poly(sp.expand(e), *X)
    return sp.Add(*[co*s**m[0]*a**m[1]*b**m[2]*c**m[3] for m, co in p.terms() if sum(m) <= k])

NS = 9
gXX = sp.series((1+s)**sp.Rational(2, 3), s, 0, NS).removeO()
gZZ = sp.series((1+s)**sp.Rational(-4, 3), s, 0, NS).removeO()
gi = [sp.Integer(-1), gXX, gZZ, gZZ]          # inverse metric, diagonal
isq = sp.series(1/(1+s), s, 0, NS).removeO()   # 1/sqrt(-g), sqrt(-g) = y = 1+s

def grad2(f, h):
    return sp.Add(*[gi[i]*sp.diff(f, X[i])*sp.diff(h, X[i]) for i in range(4)])

def box(f, k):
    return trunc(isq*sp.Add(*[sp.diff((1+s)*gi[i]*sp.diff(f, X[i]), X[i]) for i in range(4)]), k)

# ---- Synge world function: g^{ab} d_a sigma d_b sigma = 2 sigma
sig = sp.Rational(1, 2)*(-s**2 + a**2 + b**2 + c**2)
for k in range(3, DMAX+1):
    R = deg_part(trunc(grad2(sig, sig) - 2*sig, k), k)
    sig = sig - R/(2*(k-1))
assert sp.expand(trunc(grad2(sig, sig) - 2*sig, DMAX)) == 0

sa = [trunc(gi[i]*sp.diff(sig, X[i]), DMAX-1) for i in range(4)]
Bs = box(sig, 4)

# ---- U = Delta^{1/2}: 2 sigma^a U_a + (box sigma - 4) U = 0
U = sp.Integer(1)
for k in range(1, 5):
    E = trunc(2*sp.Add(*[sa[i]*sp.diff(U, X[i]) for i in range(4)]) + (Bs - 4)*U, k)
    U = U - deg_part(E, k)/(2*k)
assert sp.expand(trunc(2*sp.Add(*[sa[i]*sp.diff(U, X[i]) for i in range(4)]) + (Bs-4)*U, 4)) == 0

# ---- V0: 2 sigma^a V0_a + (box sigma - 2) V0 = - box U ;  V1(x,x) = -box V0 / 4
BU = box(U, 2)
V0 = sp.Integer(0)
for k in range(0, 3):
    E = trunc(2*sp.Add(*[sa[i]*sp.diff(V0, X[i]) for i in range(4)]) + (Bs - 2)*V0 + BU, k)
    V0 = V0 - deg_part(E, k)/(2*k+2)
V1 = sp.nsimplify(-box(V0, 0)/4)

# ---- explicit equal-time polynomials (Appendix C, Eq. hadcoef2): A = Delta_X^2, B = Delta_Z1^2 + Delta_Z2^2
Asym, Bsym = sp.symbols('A B', positive=True)
def toAB(e):
    e = sp.expand(e.subs(s, 0)).subs(c, 0); p = sp.Poly(e, a, b); out = 0
    for (i, j), co in p.terms(): out += co*Asym**(i//2)*Bsym**(j//2)
    return sp.factor(out)
print('2 sigma_4 =', toAB(2*deg_part(sig, 4)), '   2 sigma_6 =', toAB(2*deg_part(sig, 6)))
print('U_2 =', toAB(deg_part(U, 2)), '   U_4 =', toAB(deg_part(U, 4)), '   V0_2 =', toAB(deg_part(V0, 2)), '   V1(x,x) =', V1)

# ---- equal-time restriction, Laurent expansion in the spatial separation
sub0 = {s: 0}
r2h = (a**2 + b**2 + c**2)/2
sig4, sig6 = deg_part(sig, 4).subs(sub0), deg_part(sig, 6).subs(sub0)
U2, U4 = deg_part(U, 2).subs(sub0), deg_part(U, 4).subs(sub0)
V0_2 = deg_part(V0, 2).subs(sub0)
assert deg_part(V0, 0) == 0 and deg_part(V0, 1).subs(sub0) == 0
al, be = sig4/r2h, sig6/r2h
Uos_m2 = 1/r2h
Uos_0 = (U2 - al)/r2h
Uos_2 = (U4 - U2*al + al**2 - be)/r2h

pref = 1/(8*sp.pi**2)
on_sphere = {a: n, b: sp.sqrt(1 - n**2), c: 0}

# ---- adiabatic coefficients at y = 1 (massless), from the friction-free mode equation
y, u, lam = sp.symbols('y u lambda', positive=True)
gfun = sp.Function('g')(y)
om = lam*sp.sqrt(gfun)
Ssrc = 1/(4*y**2)
def Qop(W): return sp.diff(W, y, 2)/W - sp.Rational(3, 2)*(sp.diff(W, y)/W)**2
tt = sp.symbols('tt')
W2 = (Ssrc - Qop(om)/2)/(2*om)
Q1 = sp.diff(Qop(om + tt*W2), tt).subs(tt, 0)
W4 = (-Q1/2 - W2**2)/(2*om)
gy = u**2*y**sp.Rational(2, 3) + (1 - u**2)*y**sp.Rational(-4, 3)
subsd = {sp.Derivative(gfun, (y, k)): sp.diff(gy, y, k).subs(y, 1) for k in range(6, 0, -1)}
def at1(e): return sp.simplify(e.subs(subsd).subs(gfun, 1).subs(y, 1))
w2, w4 = at1(W2), at1(W4)
B2 = sp.expand(sp.simplify(-w2/(2*lam**2)*lam**3))                      # A_2 = B2/lambda^3
C4 = sp.expand(sp.simplify((-w4/(2*lam**2) + w2**2/(2*lam**3))*lam**5))  # A_4 = C4/lambda^5

def legendre_coeffs(f, lmax=10):
    return {l: sp.nsimplify(sp.Rational(2*l+1, 2)*sp.integrate(sp.expand(f*sp.legendre(l, u)), (u, -1, 1)))
            for l in range(0, lmax+1, 2)}
def I_l(l, mu):   # int_0^infty x^mu j_l(x) dx  (Abel-regularized where needed)
    return sp.sqrt(sp.pi)*sp.Integer(2)**(mu-1)*sp.gamma(sp.Rational(l+mu+1, 2))/sp.gamma(sp.Rational(l-mu+2, 2))
def FT(coeffs, mu):  # FT of sum_l c_l P_l(u) lambda^{mu-2} on the direction of separation
    return sp.Add(*[(sp.I**l)*I_l(l, mu)*cl*sp.legendre(l, n)/(2*sp.pi**2) for l, cl in coeffs.items() if cl != 0 and l > 0]) \
        + (coeffs.get(0, 0)*I_l(0, mu)/(2*sp.pi**2) if coeffs.get(0, 0) != 0 else 0)

# =========================== (i) vacuum polarization ===========================
Hfin_phi = sp.simplify((pref*Uos_0).subs(on_sphere))
cB = legendre_coeffs(B2)
assert cB[0] == 0
Afin_phi = sp.simplify(FT(cB, -1))             # FT of B2(u)/lambda^3: lambda^2 * lambda^-3 -> x^{-1}
L_phi = sp.simplify(Hfin_phi - Afin_phi)
print('B2(u) =', B2, '  Legendre:', cB)
print('Phi^2: H finite part     =', Hfin_phi)
print('Phi^2: adiabatic finite  =', Afin_phi)
print('Phi^2: H - A (should be direction independent) =', L_phi, ' ; -1/(240 pi^2) =', -1/(240*sp.pi**2))

# ===================== (ii) anisotropic combination T^X_X - T^Z_Z =====================
def op(F):  # -d_a^2 + (d_b^2 + d_c^2)/2  acting on F(Delta)  <-> <phi_X phi_X> - <phi_Z phi_Z>
    return -sp.diff(F, a, 2) + (sp.diff(F, b, 2) + sp.diff(F, c, 2))/2
Lg = (V0_2 + V1*r2h)*sp.log(r2h)
GH0 = op(Uos_2) + op(Lg)
logcoef = sp.simplify(op(V0_2 + V1*r2h))
print('V1(x,x) =', V1, '   V0_2(s=0) =', sp.factor(V0_2))
print('log coefficient in G_H (must vanish):', logcoef)
GH_const = sp.simplify((pref*GH0).subs(on_sphere))
f4 = sp.expand(sp.Rational(1, 2)*(3*u**2 - 1)*C4)
c4 = legendre_coeffs(f4)
print('C4(u) =', C4, '\nLegendre of (3u^2-1)C4/2:', c4)
assert c4[0] == 0
GA_const = sp.simplify(FT(c4, -1))           # lambda^2*lambda^-5 = lambda^-3 -> x^{-1}
D_HA = sp.simplify(GH_const - GA_const)
print('T anisotropic: H finite =', sp.expand(GH_const))
print('T anisotropic: A finite =', sp.expand(GA_const))
print('H - A (should be direction independent) =', sp.nsimplify(D_HA), ' ; +523/(498960 pi^2) =', sp.Rational(523, 498960)/sp.pi**2)

# ---- consistency of singular terms for the stress combination
GHm2 = sp.simplify((pref*op(Uos_0)).subs(on_sphere))   # Delta^-2 coefficient (unit sphere)
gB = legendre_coeffs(sp.expand(sp.Rational(1, 2)*(3*u**2 - 1)*B2))
GAm2 = sp.simplify(FT(gB, 1))
print('Delta^-2 terms equal:', sp.simplify(GHm2 - GAm2) == 0)
GHm4 = sp.simplify((pref*op(Uos_m2)).subs(on_sphere))
GAm4 = sp.simplify(FT(legendre_coeffs(sp.Rational(1, 2)*(3*u**2-1)*sp.Rational(1, 2)), 3))
print('Delta^-4 terms equal:', sp.simplify(GHm4 - GAm4) == 0)
