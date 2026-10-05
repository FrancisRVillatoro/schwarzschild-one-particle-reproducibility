import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from num_checks import propagate
M = 0.5                         # units 2M = 1
rs = np.array([1e-1,3e-2,1e-2,3e-3,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8])
wa, wb, N = 0.5, 2.5, 48
xg, wg = leggauss(N); om = 0.5*(wb-wa)*xg + 0.5*(wb+wa); W = 0.5*(wb-wa)*wg
F = np.exp(-1.0/(1.0-((om-1.5)/1.0)**2))           # smooth bump on (0.5, 2.5)
def modes(L2):
    PH, PX, BB = [], [], []
    for w in om:
        (ph_end, px_end), (xs, ph, px) = propagate(w, L2, 0j, np.sqrt(2*w)+0j, rs_out=rs, xmin=np.log(1e-12))
        o = np.argsort(np.exp(xs)); PH.append(ph[o]); PX.append(px[o])     # ordered as sorted rs
        BB.append(2*np.sqrt(M/3)*px_end)                                     # connection coefficient B_omega
    return np.array(PH), np.array(PX), np.array(BB)
rsort = np.sort(rs)
def t_of_r(r, E):
    if E == 0: return 0.0
    g = lambda x: E/((2*M/x - 1)*np.sqrt(E**2 - 1 + 2*M/x))
    return quad(g, 0, r, epsabs=1e-15, epsrel=1e-12)[0]
def run(L2, Yval, dYval, label):
    PH, PX, BB = modes(L2)
    Bstar = np.sum(W*F*BB)*Yval                      # t_* = 0
    print(f'\n{label}:  (3/2)|B*|^2 = {1.5*abs(Bstar)**2:.6e}')
    print('      r      ' + '   '.join([f'E={E:<4}' for E in (0, 0.5, 2.0)]))
    for i, r in enumerate(rsort):
        f = 2*M/r - 1; out = []
        for E in (0, 0.5, 2.0):
            t = t_of_r(r, E); ph = np.exp(-1j*om*t)
            R = np.sum(W*F*ph*PH[:, i]); Rr = np.sum(W*F*ph*PX[:, i])/r; Rt = np.sum(W*F*(-1j*om)*ph*PH[:, i])
            w, wr, wt, wth = R*Yval, Rr*Yval, Rt*Yval, R*dYval
            ur, ut = -np.sqrt(E**2 - 1 + 2*M/r), -E/f
            rho = 2*abs(ur*wr + ut*wt)**2 + (-f*abs(wr)**2 + abs(wt)**2/f + abs(wth)**2/r**2)
            out.append(r**3*rho/(1.5*abs(Bstar)**2) if abs(Bstar) > 0 else r**3*rho)
        print(f'  {r:8.0e}   ' + '   '.join([(f'{v:.6f}' if abs(Bstar)>0 else f'{v:.3e}') for v in out]))
    return Bstar
run(0, 1/np.sqrt(4*np.pi), 0.0, 'l=0 packet (Y_00)')
run(2, np.sqrt(3/(4*np.pi)), 0.0, 'l=1 packet (Y_10), endpoint at the pole theta=0')
Bp = run(2, 0.0, -np.sqrt(3/(4*np.pi)), 'l=1 packet (Y_10), endpoint on the equator (B*=0): r^3 rho shown unnormalized')
