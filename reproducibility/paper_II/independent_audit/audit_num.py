"""Numerical checks for Paper II (Annals of Physics version), units 2M=1.

Independent audit: written separately from the production scripts num_checks.py and
run_part1.py, with its own integrator (Lambert-W inversion r(eta), eta-stage from the
horizon or from Liouville-Green data, x = ln r stage down to r = 1e-12, projection on the
endpoint-Hankel pair). Numbering below refers to the Annals of Physics version of the
manuscript: Table 1, Proposition 9 (eqs. (6.3)-(6.5)), Theorem 10 (eqs. (6.8)-(6.10)),
Eq. (6.6) (directional Born benchmark). Set SKIP_A=1 to skip part (A).

 (A) Special integrals: Weber-Schafheitlin value (eq:WS), Abel value 2/pi^2, A_s and D.
 (B) Horizon-Jost modes of the massless scalar field, integrated from the horizon
     (or from first-order Liouville-Green data at K*eta = 200 for large K) down to
     r = 1e-12, projected on the endpoint-Hankel pair of the Cortez variables.
     Compared with Table 1 (s wave: kappa |beta_H|^2 ; fixed towers l=1,2,3) and
     with Theorem fixedl (alpha_H, beta_H, |beta_H|^2).

Conventions (as in the manuscript): psi = r*phi, eta = -r - ln(1-r) increasing
towards the horizon; Jost phase psi ~ exp(-i K eta + i pi/4)/sqrt(2K); endpoint
Hankel u0 = sqrt(pi*eta~/4) H_0^(2)(kappa eta~), Wronskian +i.
"""
import numpy as np, mpmath as mp
from scipy.integrate import solve_ivp, quad
from scipy.special import lambertw
import time, sys

EG = 0.5772156649015329

# ---------------- (A) special integrals ----------------
def run_A():
    mp.mp.dps = 25
    def ws_integral():
        f = lambda z: z**(-mp.mpf(1)/2)*mp.hankel2(0, z)**2
        # split [0,1] (integrable singularity) and [1,inf) oscillatory with decay z^{-3/2}
        I1 = mp.quad(f, [0, 0.01, 0.1, 1])
        I2 = mp.quadosc(f, [1, mp.inf], omega=2)   # [H0]^2 ~ exp(-2iz)
        return I1 + I2
    WS = ws_integral()
    WSclaim = -mp.sqrt(2*mp.pi)*mp.gamma(mp.mpf(1)/4)/mp.gamma(mp.mpf(3)/4)**3*(1-1j)
    print("(A) int_0^inf z^{-1/2} [H0^(2)(z)]^2 dz =", mp.nstr(WS, 12), "  claimed:", mp.nstr(WSclaim, 12))
    D = mp.pi**(mp.mpf(3)/2)*mp.gamma(mp.mpf(1)/4)/(24*mp.gamma(mp.mpf(3)/4)**3)
    As = 2*D**2
    print("    D =", mp.nstr(D, 12), " A_s = 2D^2 =", mp.nstr(As, 14), " claimed 0.417942472394")
    # beta_H leading from the first Born term: -i * (1/(6 sqrt2)) * (pi/4) * kappa^{-1/2} * WS
    betalead = -1j*(1/(6*mp.sqrt(2)))*(mp.pi/4)*WS
    print("    beta_H*sqrt(kappa) from Born term =", mp.nstr(betalead, 10), " claimed D(1+i) =", mp.nstr(D*(1+1j), 10))
    # Abel value of int z [H0^(2)]^2 dz
    def abel(eps):
        f = lambda z: z*mp.exp(-eps*z)*mp.hankel2(0, z)**2
        return mp.quad(f, [0, 0.1, 1, 5, 20, 100, 400, mp.inf])
    vals = [abel(mp.mpf(e)) for e in ('0.08', '0.04', '0.02')]
    print("    Abel int z e^{-eps z}[H0]^2: eps=0.08,0.04,0.02 ->", [mp.nstr(v, 8) for v in vals],
          " Richardson ->", mp.nstr(2*vals[2]-vals[1], 8), " claimed 2/pi^2 =", mp.nstr(2/mp.pi**2, 8))


# ---------------- (B) horizon-Jost modes ----------------
def r_of_eta(eta):
    return 1.0 + lambertw(-np.exp(-1.0 - eta)).real

def Vpot_eta(r, L2):
    """potential in psi_eta_eta + [K^2 + V] psi = 0"""
    return (1 - r)/r**4 + L2*(1 - r)/r**3

def dV_deta(r, L2):
    dVdr = -(1/r**4) - 4*(1 - r)/r**5 + L2*(-(1/r**3) - 3*(1 - r)/r**4)
    return dVdr*(1 - r)/r

def jost_data_eta(K, L2, r0):
    """first-order Liouville-Green Jost data (psi, psi_eta) at r0, with the exact
    asymptotic phase: psi -> exp(-i K eta + i pi/4)/sqrt(2K) as eta -> inf."""
    eta0 = -r0 - np.log(1 - r0)
    W0 = np.sqrt(K**2 + Vpot_eta(r0, L2))
    # I0 = int_{eta0}^{inf} (W-K) d eta = int_{r0}^{1} V/(W+K) * r/(1-r) dr
    g = lambda r: Vpot_eta(r, L2)/(np.sqrt(K**2 + Vpot_eta(r, L2)) + K)*r/(1 - r)
    I0, err = quad(g, r0, 1.0, limit=400, epsabs=1e-14, epsrel=1e-13)
    psi0 = (2*W0)**-0.5*np.exp(-1j*K*eta0 + 1j*I0 + 1j*np.pi/4)
    dpsi0 = (-1j*W0 - dV_deta(r0, L2)/(4*W0**2))*psi0
    return eta0, psi0, dpsi0, I0

def integrate_jost(K, L2, r_start, r_switch=0.3, r_min=1e-12, rtol=1e-12, atol=1e-16):
    """returns (a, b) with phi = a + b ln r + ... at r_min, plus diagnostics"""
    eta0, psi0, dpsi0, I0 = jost_data_eta(K, L2, r_start)
    if r_start > r_switch:
        # stage 1 in eta from eta0 down to eta(r_switch)
        def rhs_eta(eta, y):
            r = r_of_eta(eta)
            return [y[1], -(K**2 + Vpot_eta(r, L2))*y[0]]
        eta1 = -r_switch - np.log(1 - r_switch)
        sol = solve_ivp(rhs_eta, [eta0, eta1], [psi0, dpsi0], method='DOP853', rtol=rtol, atol=atol)
        psi1, dpsi1 = sol.y[0, -1], sol.y[1, -1]
        r1 = r_switch
    else:
        psi1, dpsi1, r1 = psi0, dpsi0, r_start
    # convert to phi, phi_x at r1
    phi1 = psi1/r1
    dphidx1 = dpsi1*r1/(1 - r1) - psi1/r1      # psi_r = psi_eta * r/(1-r); phi_x = psi_r - psi/r
    def rhs_x(x, y):
        r = np.exp(x)
        Vr = L2 + K**2*r**3/(1 - r)
        return [y[1], (r*y[1] - r*Vr*y[0])/(1 - r)]
    sol2 = solve_ivp(rhs_x, [np.log(r1), np.log(r_min)], [phi1, dphidx1], method='DOP853', rtol=rtol, atol=atol)
    phi, dphi = sol2.y[0, -1], sol2.y[1, -1]
    b = dphi
    a = phi - b*np.log(r_min)
    # Wronskian check: r^2 (phi phi_eta* - c.c.) = i  with phi_eta = phi_x * (1-r)/r
    r = r_min
    wr = (1 - r)*(phi*np.conj(dphi) - dphi*np.conj(phi))   # = psi psi_eta^* - psi_eta psi^* = i for a normalized mode
    return a, b, wr, I0

def endpoint_coeffs(a, b, K, L2):
    """Bogoliubov coefficients of the solution (a,b) relative to the endpoint-Hankel pair."""
    L = np.sqrt(L2)
    kappa = np.sqrt(K**2 + L2)
    if L2 == 0:
        A = np.sqrt(2)*(a + 0.5*b*np.log(2))
        B = b/np.sqrt(2)
    else:
        lhat = L/kappa
        B = np.sqrt(2)*b
        A = (a - 2*b*np.log(2*lhat))/np.sqrt(2)
    ar = np.sqrt(np.pi)/2
    ai = -(np.log(kappa/2) + EG)/np.sqrt(np.pi)
    amb = 1j*np.sqrt(np.pi)*B
    apb = (A - 1j*ai*amb)/ar
    alpha = (apb + amb)/2
    beta = (apb - amb)/2
    return alpha, beta, kappa

def formula(K, l):
    Lam = K/(4*np.exp(3*EG)*(l*(l+1))**2)
    al = 1.25 - 1j/(2*np.pi)*np.log(Lam)
    be = -0.75 - 1j/(2*np.pi)*np.log(Lam)
    return al, be, 9/16 + np.log(Lam)**2/(4*np.pi**2), -np.log(Lam)/(2*np.pi)

if __name__ == '__main__':
    import os
    if not os.environ.get('SKIP_A'):
        run_A()
    print("\n(B1) massless s wave: kappa |beta_H|^2 (Table 1, top); start at r0 = 1-1e-9 (horizon side)")
    paper = {2: 0.306176, 4: 0.333939, 8: 0.356089, 16: 0.372997, 32: 0.385566, 64: 0.394755, 128: 0.401401}
    for kap in [2, 4, 8, 16, 32, 64, 128]:
        t0 = time.time()
        a, b, wr, I0 = integrate_jost(kap, 0.0, 1 - 1e-9)
        al, be, kk = endpoint_coeffs(a, b, kap, 0.0)
        print(f"  kappa={kap:4d}: kappa|beta|^2 = {kap*abs(be)**2:.6f}  (paper {paper[kap]:.6f}); |alpha|^2-|beta|^2-1 = {abs(al)**2-abs(be)**2-1:.1e}; "
              f"Wronskian = {wr:.6f}; beta*sqrt(kappa) = {be*np.sqrt(kap):.4f}; {time.time()-t0:.1f}s")
        sys.stdout.flush()
    a, b, wr, I0 = integrate_jost(512, 0.0, 1 - 1e-9)
    al, be, kk = endpoint_coeffs(a, b, 512, 0.0)
    print(f"  kappa=512: beta*sqrt(kappa) = {be*np.sqrt(512):.4f} (paper 0.4481+0.4569i), kappa|beta|^2 = {512*abs(be)**2:.6f}")

    print("\n(B2) validation of the Liouville-Green start at K*eta=200 (l=1, K=1e4): horizon start vs LG start")
    K = 1e4; L2 = 2.0
    a1, b1, wr1, _ = integrate_jost(K, L2, 1 - 1e-9)
    al1, be1, _ = endpoint_coeffs(a1, b1, K, L2)
    eta0 = 200/K; r0 = r_of_eta(eta0)
    a2, b2, wr2, _ = integrate_jost(K, L2, r0)
    al2, be2, _ = endpoint_coeffs(a2, b2, K, L2)
    print(f"  horizon start: alpha={al1:.5f} beta={be1:.5f} |beta|^2={abs(be1)**2:.6f}")
    print(f"  LG start r0={r0:.4f}: alpha={al2:.5f} beta={be2:.5f} |beta|^2={abs(be2)**2:.6f}")

    print("\n(B3) fixed angular towers (Table 1, bottom) with LG start at K*eta=200")
    rows = [(1, 1e4, 1.00864, 1.2489-0.6700j, -0.6871-0.7325j), (1, 1e5, 1.75205, 1.2536-1.0866j, -0.7262-1.1067j),
            (1, 1e6, 2.73394, 1.2523-1.4716j, -0.7413-1.4780j), (1, 1e7, 3.96843, 1.2511-1.8448j, -0.7469-1.8468j),
            (2, 1e5, 1.04792, 1.2497-0.6973j, -0.6954-0.7512j), (2, 1e7, 2.78708, 1.2520-1.4898j, -0.7425-1.4953j),
            (3, 1e5, 0.70965, 1.2366-0.4248j, -0.6565-0.5279j), (3, 1e7, 2.16185, 1.2528-1.2619j, -0.7365-1.2726j)]
    for l, K, nb_paper, al_paper, be_paper in rows:
        L2 = l*(l+1)
        r0 = r_of_eta(200/K)
        t0 = time.time()
        a, b, wr, _ = integrate_jost(K, L2, r0)
        al, be, kk = endpoint_coeffs(a, b, K, L2)
        alf, bef, nbf, imf = formula(K, l)
        print(f"  l={l} K={K:.0e}: |beta|^2={abs(be)**2:.5f} (paper {nb_paper:.5f}, formula {nbf:.5f}); alpha={al:.4f} (paper {al_paper:.4f}, formula {alf:.4f}); "
              f"beta={be:.4f} (paper {be_paper:.4f}, formula {bef:.4f}); norm err {abs(al)**2-abs(be)**2-1:.1e}; {time.time()-t0:.1f}s")
        sys.stdout.flush()

    print("\n(B4) directional Born benchmark along rays lhat^2 = 1/2 and 4/5 (Eq. (6.6), |beta_H|^2 only)")
    for lhat2 in (0.5, 0.8):
        for kap in (10.0, 40.0):
            L2 = lhat2*kap**2
            K = np.sqrt(kap**2 - L2)
            a, b, wr, _ = integrate_jost(K, L2, 1 - 1e-9)
            al, be, kk = endpoint_coeffs(a, b, K, L2)
            bench = 1/(36*np.pi**2*lhat2**2)
            val = abs(be)**2*kap**4
            print(f"  lhat^2={lhat2}, kappa={kap}: |beta|^2 kappa^4 = {val:.5f}, benchmark 1/(36 pi^2 lhat^4) = {bench:.5f}, rel dev {100*(val/bench-1):+.2f}%")
            sys.stdout.flush()
