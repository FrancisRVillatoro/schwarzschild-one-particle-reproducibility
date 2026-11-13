"""Phase 0 for a fermionic paper III (units 2M = 1). Conventions of the manuscript: H = gamma^0[gamma^1 k/a + gamma-tilde^2 lambda/r + m].
Separated Dirac mode system in the Schwarzschild interior (Kantowski-Sachs), in the Dirac-algebra
representation unitarily equivalent to Chandrasekhar's separated system:
    i d chi/d sigma = H(sigma) chi,   H = gamma^0 [ gamma^1 (k/a) + gamma^2 (kappa/r) + m ],
with a = sqrt(f), f = 1/r - 1, kappa = +-(j+1/2), k the Killing frequency (integer in the compactified theory),
chi = (a r^2)^{1/2} psi (volume rescaling). H is Hermitian, so the norm of chi is conserved.
Checks: (1) chi has a limit at r -> 0 (no logarithm), rate; (2) the instantaneous (Parker) occupation converges;
(3) endpoint Bogoliubov block relative to the transverse projector: |beta|^2 as a function of (k, kappa) and the
scaling variable kappa / k^{1/4} that controls adiabatic following of the spin rotation."""
import numpy as np
from scipy.integrate import solve_ivp
I2=np.eye(2); s1=np.array([[0,1],[1,0]],complex); s2=np.array([[0,-1j],[1j,0]]); s3=np.diag([1.,-1.]).astype(complex)
g0=np.block([[I2,0*I2],[0*I2,-I2]]).astype(complex)
def gam(s): return np.block([[0*I2,s],[-s,0*I2]])
g1=gam(s1); g2=np.block([[0*I2,-s3],[s3,0*I2]])   # gamma-tilde^2 of the manuscript, (gamma^0,gamma^1,gamma-tilde^2) is a Clifford triple
def H(r,k,kap,m):
    a=np.sqrt(1/r-1)
    return g0@(g1*(k/a)+g2*(kap/r)+m*np.eye(4))
def proj(Hm,sign):
    w,v=np.linalg.eigh(Hm); sel=v[:,w*sign>0]; return sel@sel.conj().T
def evolve(k,kap,m,r0,rs):
    """propagate the two positive-energy eigenvectors of H(r0); return chi(r) at the radii rs (sorted descending)."""
    H0=H(r0,k,kap,m); w,v=np.linalg.eigh(H0); V0=v[:,w>0]            # 4x2
    def rhs(r,y):
        Y=(y[:8]+1j*y[8:]).reshape(4,2); dY=-1j*(H(r,k,kap,m)/np.sqrt(1/r-1))@Y
        return np.concatenate([dY.real.ravel(),dY.imag.ravel()])
    y0=np.concatenate([V0.real.ravel(),V0.imag.ravel()])
    sol=solve_ivp(rhs,(r0,rs[-1]),y0,t_eval=rs,method='DOP853',rtol=1e-11,atol=1e-13)
    return [ (sol.y[:8,i]+1j*sol.y[8:,i]).reshape(4,2) for i in range(len(rs)) ]
if __name__ == '__main__':
    rs=np.array([0.5,0.2,0.1,0.03,0.01,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8])
    print('=== (1)-(2) convergence at the singularity, k=5, kappa=2, m=0, start r0=0.9 ===')
    Ys=evolve(5,2,0.,0.9,rs); Yend=Ys[-1]
    Pm_end=proj(g0@g2*np.sign(2),-1)      # transverse endpoint projector onto negative energy
    for r,Y in zip(rs,Ys):
        n_parker=np.trace(Y.conj().T@proj(H(r,5,2,0.),-1)@Y).real/2   # occupation per initial positive state
        print(f'  r={r:8.0e}  norm-1={np.linalg.norm(Y)**2/2-1:+.1e}  |chi(r)-chi(end)|={np.linalg.norm(Y-Yend):.3e}  ratio to r^1/2={np.linalg.norm(Y-Yend)/r**0.5:.3f}  Parker n={n_parker:.6f}  endpoint beta^2={np.trace(Y.conj().T@Pm_end@Y).real/2:.6f}')
    print('\n=== (3) endpoint Bogoliubov block vs (k, kappa): |beta|^2 = tr(P_+(r0) U^+ P_-^perp U P_+(r0)) / 2 ===')
    ks=[3,10,30,100,300,1000]; kaps=[1,2,3,5,8,12,20,30,50,80]
    rend=np.array([1e-6])
    tab={}
    for k in ks:
        row=[]
        for kap in kaps:
            Y=evolve(k,kap,0.,0.9,rend)[0]; b2=np.trace(Y.conj().T@proj(g0@g2,-1)@Y).real/2; row.append(b2); tab[(k,kap)]=b2
        print('k=%5d  '%k+'  '.join(f'{b:.4f}' for b in row))
    print('kappa:        '+'  '.join(f'{kap:6d}' for kap in kaps))
    print('\nscaling variable x = kappa / k^{1/4} and |beta|^2:')
    pts=sorted([(kap/k**0.25,b,k,kap) for (k,kap),b in tab.items()])
    for x,b,k,kap in pts:
        if b>1e-3: print(f'  x={x:6.3f}  |beta|^2={b:.4f}   (k={k}, kappa={kap})')
    print('\nsensitivity to the starting slice r0 (k=100, kappa=3): ')
    for r0 in (0.6,0.9,0.99):
        Y=evolve(100,3,0.,r0,rend)[0]; print(f'  r0={r0}: |beta|^2={np.trace(Y.conj().T@proj(g0@g2,-1)@Y).real/2:.4f}')
    print('massive case m=1 (k=30, kappa=3):', end=' ')
    Y=evolve(30,3,1.,0.9,rend)[0]; print(f'|beta|^2={np.trace(Y.conj().T@proj(g0@g2,-1)@Y).real/2:.4f}  (massless {tab[(30,3)]:.4f})')
