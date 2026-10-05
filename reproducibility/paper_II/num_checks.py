"""
Numerical checks for paper II (units 2M = 1, so M = 1/2).

Part 1. Fixed angular towers: exact horizon-Jost solution of
    psi_etaeta + [K^2 + (1-r)(1 + L^2 r)/r^4] psi = 0,   eta = -r - ln(1-r),  psi = r phi,
with psi ~ exp(-i(K eta - pi/4))/sqrt(2K) at the horizon, propagated to r -> 0 and projected
onto the endpoint-Hankel pair; compared with Theorem (fixed-l obstruction):
    alpha_H = 5/4 - (i/2pi) ln Lambda,  beta_H = -3/4 - (i/2pi) ln Lambda,
    Lambda = K / (4 e^{3 gamma} [l(l+1)]^2).
The massless s wave is used as a validation against the published table.

Part 2. One-particle stress: conditioned interior packet built from right-interior modes
(psi ~ exp(-i omega r_*) at the horizon), energy density along radial geodesics of
Killing energy E, compared with the limit (3/2)|B_*|^2 / r^3.
"""
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

EG = np.euler_gamma
RT, AT = 1e-11, 1e-13

def p_of_eta(eta):
    if eta > 5:
        p = np.exp(-(eta + 1.0))
        for _ in range(60):
            p = np.exp(-(eta + 1.0 - p))
        return p
    r = brentq(lambda r: -r - np.log1p(-r) - eta, 1e-14, 1 - 1e-15, xtol=1e-15)
    return 1 - r

def eta_of_r(r):
    return -r - np.log1p(-r)

def propagate(K, L2, cm0, cp0, r1=0.3, eta_max=None, rs_out=None, xmin=np.log(1e-10)):
    """Return (phi, phi_x) at xmin (and at rs_out if given), phi_x = r dphi/dr."""
    if eta_max is None:
        eta_max = 30.0 + np.log(1 + L2)
    eta1 = eta_of_r(r1)
    def f_eta(eta, y):
        cm, cp, p = y[0] + 1j*y[1], y[2] + 1j*y[3], np.exp(y[4])
        r = 1 - p
        V = p*(1 + L2*r)/r**4
        e = np.exp(2j*K*eta)
        dcp = 1j*V/(2*K)*(cm/e + cp)
        dcm = -1j*V/(2*K)*(cm + cp*e)
        return [dcm.real, dcm.imag, dcp.real, dcp.imag, -1.0/r]
    y0 = [cm0.real, cm0.imag, cp0.real, cp0.imag, np.log(p_of_eta(eta_max))]
    s1 = solve_ivp(f_eta, (eta_max, eta1), y0, method='DOP853', rtol=RT, atol=AT)
    cm, cp, p = s1.y[0, -1] + 1j*s1.y[1, -1], s1.y[2, -1] + 1j*s1.y[3, -1], np.exp(s1.y[4, -1])
    r = 1 - p
    psi = (cm*np.exp(-1j*K*eta1) + cp*np.exp(1j*K*eta1))/np.sqrt(2*K)
    dpsi = (-1j*K*cm*np.exp(-1j*K*eta1) + 1j*K*cp*np.exp(1j*K*eta1))/np.sqrt(2*K)
    phi, phix = psi/r, r*dpsi/(1 - r) - psi/r
    def f_x(x, y):
        rr = np.exp(x)
        ph, px = y[0] + 1j*y[1], y[2] + 1j*y[3]
        pxx = (rr*px - rr*L2*ph - K**2*rr**4*ph/(1 - rr))/(1 - rr)
        return [px.real, px.imag, pxx.real, pxx.imag]
    teval = None
    if rs_out is not None:
        teval = np.sort(np.log(rs_out))[::-1]
    s2 = solve_ivp(f_x, (np.log(r), xmin), [phi.real, phi.imag, phix.real, phix.imag],
                   method='DOP853', rtol=RT, atol=AT, t_eval=teval)
    out_end = (s2.y[0, -1] + 1j*s2.y[1, -1], s2.y[2, -1] + 1j*s2.y[3, -1]) if teval is None else None
    if teval is not None:
        s3 = solve_ivp(f_x, (np.log(r), xmin), [phi.real, phi.imag, phix.real, phix.imag],
                       method='DOP853', rtol=RT, atol=AT)
        out_end = (s3.y[0, -1] + 1j*s3.y[1, -1], s3.y[2, -1] + 1j*s3.y[3, -1])
        prof = (s2.t, s2.y[0] + 1j*s2.y[1], s2.y[2] + 1j*s2.y[3])
        return out_end, prof
    return out_end

def endpoint_coeffs(K, l, xmin=np.log(1e-10)):
    L2 = l*(l + 1)
    (ph, px) = propagate(K, L2, np.exp(1j*np.pi/4), 0j, xmin=xmin)
    b = px
    a = ph - b*xmin
    if l == 0:
        av = np.sqrt(np.pi/8)*(1 - 2j/np.pi*(np.log(K/4) + EG))
        bv = -1j*np.sqrt(2/np.pi)
    else:
        L = np.sqrt(L2)
        av = np.sqrt(np.pi/2)*(1 - 2j/np.pi*(np.log(L) + EG))
        bv = -1j/np.sqrt(2*np.pi)
    A = np.array([[av, np.conj(av)], [bv, np.conj(bv)]])
    al, be = np.linalg.solve(A, np.array([a, b]))
    return al, be

# ---------------- large-K start: first-order Liouville-Green data at r1 (K*eta1 = 200) ----------------
from scipy.integrate import quad
def endpoint_coeffs_LG(K, l, xmin=np.log(1e-12)):
    L2 = l*(l + 1)
    r1 = min(0.3, np.sqrt(400.0/K))
    V = lambda r: (1 - r)*(1 + L2*r)/r**4
    I = quad(lambda r: (1 + L2*r)/(r**3*(np.sqrt(K**2 + V(r)) + K)), r1, 1, limit=500, epsabs=1e-15, epsrel=1e-13)[0]
    eta1, Q = eta_of_r(r1), K**2 + V(r1)
    h = 1e-6*r1; dVdeta = (V(r1 + h) - V(r1 - h))/(2*h)*(1 - r1)/r1
    psi = np.exp(1j*np.pi/4)*(2*np.sqrt(Q))**-0.5*np.exp(-1j*(K*eta1 - I))
    dpsi = (-1j*np.sqrt(Q) - dVdeta/(4*Q))*psi
    phi, phix = psi/r1, r1*dpsi/(1 - r1) - psi/r1
    def f_x(x, y):
        rr = np.exp(x); ph, px = y[0] + 1j*y[1], y[2] + 1j*y[3]
        pxx = (rr*px - rr*L2*ph - K**2*rr**4*ph/(1 - rr))/(1 - rr)
        return [px.real, px.imag, pxx.real, pxx.imag]
    sx = solve_ivp(f_x, (np.log(r1), xmin), [phi.real, phi.imag, phix.real, phix.imag], method='DOP853', rtol=RT, atol=AT)
    ph, px = sx.y[0, -1] + 1j*sx.y[1, -1], sx.y[2, -1] + 1j*sx.y[3, -1]
    b = px; a = ph - b*xmin
    if l == 0:
        av = np.sqrt(np.pi/8)*(1 - 2j/np.pi*(np.log(K/4) + EG)); bv = -1j*np.sqrt(2/np.pi)
    else:
        L = np.sqrt(L2); av = np.sqrt(np.pi/2)*(1 - 2j/np.pi*(np.log(L) + EG)); bv = -1j/np.sqrt(2*np.pi)
    return np.linalg.solve(np.array([[av, np.conj(av)], [bv, np.conj(bv)]]), np.array([a, b]))

if __name__ == '__main__':
    print('=== Massless s wave: kappa |beta_H|^2 (Table 1, top) ===')
    for K, ref in [(2, 0.306176), (4, 0.333939), (8, 0.356088), (16, 0.372996), (32, 0.385565), (64, 0.394754), (128, 0.401400)]:
        al, be = endpoint_coeffs(K, 0)
        print(f'kappa={K:4d}  kappa|beta|^2={K*abs(be)**2:.6f}  (published {ref:.6f})  |alpha|^2-|beta|^2-1={abs(al)**2-abs(be)**2-1:+.1e}')
    print('Fixed angular towers: run_part1.py.  Packet convergence: run_part2.py.')
