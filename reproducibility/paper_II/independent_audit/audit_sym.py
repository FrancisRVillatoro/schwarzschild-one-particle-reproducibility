"""Symbolic checks for Paper II (Annals of Physics version).

Independent audit, written separately from paper2_symbolic_checks.py. Numbering below
refers to the Annals of Physics version of the manuscript: Eq. (2.3) (pump), Lemma 1,
Eq. (3.2), (3.4), (3.6) (count prefactors), Eq. (3.5) (Cortez variables), Lemma 5 and
Eq. (6.2) (deviations and s-wave series), Eq. (8.9) (radial geodesics), Appendix F
(Frobenius coefficient), Remark 18 (curvature coupling).

Units 2M=1. Checks:
 (1) sigma(r), q(r), H(r), the pump Q(r) of Eq. (eq:Q), and Q = H'/2 + H^2/4.
 (2) Cortez-type rescaling b_lhat, harmonic time, conformal time eta-tilde,
     mass function s_lhat (Appendix A formula), and the deviations V_1
     of Lemma dV: l=0 (s_0=(1-r)/r^4, V_1=(1/3)r^-3[1+O(r)], eta-tilde series,
     the s-wave expansion (eq:Vs)), and l>=1 (V_1 -> 1/(3 lhat^2)).
 (3) lambda (b' = lambda b) for l=0 and l>=1.
 (4) Radial geodesics: sigma_E(r) and t(r) expansions (eq:geodesic), with M restored.
 (5) Frobenius coefficient d = -l(l+1) b/(2M) for phi = a + b ln r + r(c + d ln r).
 (6) Parker count prefactor and velocity count prefactor 9/8.
"""
import sympy as sp

r = sp.symbols('r', positive=True)
print("== (1) geometry, 2M=1")
f = 1/r - 1
sigma = sp.asin(sp.sqrt(r)) - sp.sqrt(r*(1-r))
dsig = sp.simplify(sp.diff(sigma, r) - 1/sp.sqrt(f))
print("d sigma/dr - 1/sqrt(f) =", dsig)
print("sigma series:", sp.series(sigma, r, 0, 4))
q = r**2*sp.sqrt(f)
print("q =", sp.simplify(q), " series:", sp.series(sp.simplify(q), r, 0, 3))
# derivative wrt sigma: d/dsigma = sqrt(f) d/dr
D = lambda expr: sp.sqrt(f)*sp.diff(expr, r)
H = sp.simplify(D(q)/q)
print("H =", sp.simplify(H), "  claimed (3-4r)/(2 r^{3/2} sqrt(1-r)):",
      sp.simplify(H - (3-4*r)/(2*r**sp.Rational(3,2)*sp.sqrt(1-r))))
W = sp.sqrt(q)
Qpump = sp.simplify(D(D(W))/W)
Qclaim = -(1/r**3)*(sp.Rational(1,2) + 1/(16*(1-r)))
print("Q - claimed =", sp.simplify(Qpump - Qclaim))
print("Q - (H'/2 + H^2/4) =", sp.simplify(Qpump - (D(H)/2 + H**2/4)))
# Q in terms of sigma near 0: Q*sigma^2 -> -1/4
print("Q*sigma^2 series:", sp.series(sp.simplify(Qpump*sigma**2), r, 0, 2))

print("\n== (2) Cortez rescaling, l=0")
# b_0 = r^2, w = r(1-r), harmonic time dtau_h = dsigma/q ; eta-tilde = int b/w dr
b0 = r**2
w = r*(1-r)
L0 = sp.diff(b0, r)/b0
s0 = sp.simplify(w*(w*L0**2 - 2*sp.diff(w, r)*L0 - 2*w*sp.diff(L0, r))/(4*b0**2))
print("s_0 =", s0, " claimed (1-r)/r^4:", sp.simplify(s0 - (1-r)/r**4))
eta0 = sp.integrate(b0/w, (r, 0, r))
print("eta-tilde(l=0) =", sp.simplify(eta0), " series:", sp.series(-r - sp.log(1-r), r, 0, 6))
V1 = s0 - 1/(4*(-r - sp.log(1-r))**2)
print("V_1 (l=0) series in r:", sp.series(V1, r, 0, 1))
# express V_1 in powers of eta-tilde: invert eta = r^2/2 + r^3/3 + r^4/4 + ...
eta = sp.symbols('eta', positive=True)
# r as series in sqrt(2 eta): r = x - x^2/3 + ... solve iteratively
x = sp.sqrt(2*eta)
rs = x
for _ in range(6):
    rs = sp.series(sp.sqrt(2*(eta - (sp.series(-rs - sp.log(1-rs), eta, 0, 4).removeO() - rs**2/2))), eta, 0, 3).removeO()
rs = sp.simplify(rs)
print("r(eta) =", sp.series(rs, eta, 0, 3))
V1eta = sp.series(V1.subs(r, rs), eta, 0, 1)
print("V_1 as series in eta-tilde:", V1eta)
print("claimed coefficients: 1/(6 sqrt2) =", sp.nsimplify(1/(6*sp.sqrt(2))), "=", float(1/(6*sp.sqrt(2))),
      "; -23/(540 sqrt2) =", float(-23/(540*sp.sqrt(2))))
# numeric check of the claimed expansion directly
import mpmath as mp
mp.mp.dps = 30
def V1_num(rv):
    rv = mp.mpf(rv)
    s0v = (1-rv)/rv**4
    etv = -rv - mp.log(1-rv)
    return s0v - 1/(4*etv**2), etv
for rv in ['1e-2','1e-3','1e-4']:
    v, et = V1_num(rv)
    lead = 1/(6*mp.sqrt(2))*et**(-mp.mpf(3)/2)
    sec = -mp.mpf(23)/(540*mp.sqrt(2))*et**(-mp.mpf(1)/2)
    print(f"r={rv}: V1={mp.nstr(v,12)}, lead={mp.nstr(lead,12)}, lead+sec={mp.nstr(lead+sec,12)}, resid/(1)={mp.nstr(v-lead-sec,6)}")

print("\n== (2b) l>=1 rescaling")
lh = sp.symbols('lhat', positive=True)
b = sp.sqrt((1-lh**2)*r**4 + lh**2*r*(1-r))
Lb = sp.diff(b, r)/b
s_l = w*(w*Lb**2 - 2*sp.diff(w, r)*Lb - 2*w*sp.diff(Lb, r))/(4*b**2)
print("s_lhat series:", sp.series(sp.simplify(s_l), r, 0, 1))
print("claimed (16 lhat^2 r)^{-1}[1+5r]: coefficient check ->",
      sp.simplify(sp.series(sp.simplify(s_l)*16*lh**2*r, r, 0, 2)))
# eta-tilde for l>=1: integral of b/w from 0 to r, series
integrand = sp.series(b/w, r, 0, 3).removeO()
etal = sp.integrate(integrand, (r, 0, r))
print("eta-tilde (l>=1) series:", sp.simplify(etal), " = 2 lhat sqrt(r)[1 + r/6 + ...]?",
      sp.series(sp.simplify(etal/(2*lh*sp.sqrt(r))), r, 0, 2))
# V1 = s - 1/(4 eta^2)
etal_full = sp.integrate(sp.series(b/w, r, 0, 5).removeO(), (r, 0, r))
V1l = sp.series(sp.simplify(s_l) - 1/(4*etal_full**2), r, 0, 1)
print("V_1 (l>=1) series:", sp.simplify(V1l))
print("claimed 1/(3 lhat^2) ->", sp.simplify(sp.limit(sp.simplify(s_l) - 1/(4*etal_full**2), r, 0) - 1/(3*lh**2)))

print("\n== (3) lambda: b' = lambda b with d/dtau_h = q d/dsigma = q sqrt(f) d/dr = w d/dr")
for name, bb in [('l=0', b0), ('l>=1', b)]:
    lam = sp.series(sp.simplify(w*sp.diff(bb, r)/bb), r, 0, 1)
    print(name, " b'/b =", lam)

print("\n== (4) geodesics with M restored")
M, E = sp.symbols('M E', positive=True)
rr = sp.symbols('r', positive=True)
sigE = sp.integrate(sp.series((2*M/rr + E**2 - 1)**(-sp.Rational(1,2)), rr, 0, 3).removeO(), (rr, 0, rr))
print("sigma_E =", sp.simplify(sigE))
print("claimed 2 r^{3/2}/(3 sqrt(2M)) [1 - 3(E^2-1) r/(20M)]: ratio series ->",
      sp.series(sp.simplify(sigE/(2*rr**sp.Rational(3,2)/(3*sp.sqrt(2*M)))), rr, 0, 2))
dtdr = E/((2*M/rr - 1)*sp.sqrt(2*M/rr + E**2 - 1))
tser = sp.integrate(sp.series(dtdr, rr, 0, 4).removeO(), (rr, 0, rr))
lead_t = sp.limit(tser/rr**sp.Rational(5,2), rr, 0)
print("t(r)-t_* =", sp.simplify(tser), "  leading coefficient:", sp.simplify(lead_t),
      "  claimed 2E/(5 (2M)^{3/2}) ->", sp.simplify(lead_t - 2*E/(5*(2*M)**sp.Rational(3,2))))
print("sigma_E^2 leading = 2 r^3/(9M)?", sp.simplify(sp.limit(sigE**2/(2*rr**3/(9*M)), rr, 0)))
# v_E
print("v_E = -E (2M/r - 1 + E^2)^{-1/2} ~ -E (r/2M)^{1/2}:", sp.series(-E*(2*M/rr - 1 + E**2)**(-sp.Rational(1,2))/( -E*sp.sqrt(rr/(2*M))), rr, 0, 2))

print("\n== (5) Frobenius: d = -l(l+1) b/(2M)")
l_, om, a_, b_, c_, d_ = sp.symbols('l omega a b c d')
phi = a_ + b_*sp.log(rr) + rr*(c_ + d_*sp.log(rr))
expr = sp.diff(rr*(2*M - rr)*sp.diff(phi, rr), rr) + (l_*(l_+1) + om**2*rr**3/(2*M - rr))*phi
expr = sp.expand(expr)
# collect terms of order r^0 (times log and const)
expr0 = sp.series(expr, rr, 0, 1).removeO()
coef_log = sp.simplify(expr0.coeff(sp.log(rr)))
coef_const = sp.simplify(expr0 - coef_log*sp.log(rr))
print("order r^0, log r coefficient:", coef_log, "-> d =", sp.solve(coef_log, d_))
print("order r^0, constant:", coef_const, "-> c =", sp.solve(coef_const, c_))

print("\n== (6) counts")
# Parker: n+1/2 = (|p|^2/(2q) + q w^2|phi|^2/2)/w ; with p = sqrt(3/2) B, q = 3 sigma/2 : |p|^2/(2q w) = |B|^2/(2 w sigma)
print("Parker prefactor: |B|^2/(2 omega sigma)  [ (3/2)|B|^2 / (2*(3/2)sigma*omega) ] ok")
# velocity: n_chi ~ (|B|^2/2) ln^2 sigma ; ln sigma = (3/2) ln r + const -> (9/8)|B|^2 ln^2 r
print("velocity prefactor: (1/2)*(3/2)^2 =", sp.Rational(1,2)*sp.Rational(9,4))
# n_ref: lambda^2 |B|^2 ln^2 sigma /(12 kappa b): from |Pi~|^2/(2 kappa) with Pi~ = lambda phi/(2 sqrt b), phi = sqrt(2/3) B ln sigma
lam, B, kap, bb = sp.symbols('lambda B kappa b', positive=True)
Pit = lam*sp.sqrt(sp.Rational(2,3))*B/(2*sp.sqrt(bb))
print("n_ref prefactor:", sp.simplify(Pit**2/(2*kap)), " claimed lambda^2 B^2/(12 kappa b)")

print("\n== (7) improvement term on Kasner: Phi^2 = C ln^2 sigma, volume ~ sigma")
s_ = sp.symbols('sigma', positive=True); C = sp.symbols('C', positive=True)
Phi2 = C*sp.log(s_)**2
box = -(1/s_)*sp.diff(s_*sp.diff(Phi2, s_), s_)
rho_imp = sp.simplify(-box - sp.diff(Phi2, s_, 2))
print("rho_imp/xi =", rho_imp, "  claimed 2C ln(sigma)/sigma^2")
