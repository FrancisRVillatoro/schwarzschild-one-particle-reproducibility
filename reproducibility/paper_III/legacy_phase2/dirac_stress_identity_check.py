"""Check of the exact energy-density identity of Theorem F (paper III), with random data (numpy, 20 draws).
Psi = sum_k c_k e^{ikt} Psi_k,  Psi_k = (F1 Theta + F2 s1 Theta, G1 Theta + G2 s1 Theta),  i d/dsigma chi_k = H_k chi_k.
Identity:  Re[ Psi^+ (i d_sigma Psi) ]  =  sum_{k,k'} Re[ conj(c_k) c_k' e^{i(k'-k)t} ( chi_k^+ Hbar_{kk'} chi_k' |Theta|^2 + Xi_{kk'} Theta^+ s1 Theta ) ],
   Hbar_{kk'} = H at (k+k')/2,   Xi_{kk'} = (k+k')/(2a) chi_k^+ g5 chi_k' + m chi_k^+ g0 Sigma1 chi_k',
with g5 = [[0,I],[I,0]], Sigma1 = diag(s1,s1).  Theta is an arbitrary 2-spinor and chi_k arbitrary 4-vectors (no equation needed)."""
import numpy as np
rng=np.random.default_rng(7)
s1=np.array([[0,1],[1,0]],complex); s3=np.diag([1.,-1.]).astype(complex); I2=np.eye(2); Z=np.zeros((2,2))
g0=np.block([[I2,Z],[Z,-I2]]).astype(complex); g1=np.block([[Z,s1],[-s1,Z]]); g2t=np.block([[Z,-s3],[s3,Z]])
g5=np.block([[Z,I2],[I2,Z]]).astype(complex); Sig1=np.block([[s1,Z],[Z,s1]])
def H(k,a,lam,r,m): return g0@(g1*k/a+g2t*lam/r+m*np.eye(4))
def Psi(chi,th): return np.concatenate([chi[0]*th+chi[1]*(s1@th), chi[2]*th+chi[3]*(s1@th)])
worst=0.
for trial in range(20):
    ks=rng.normal(size=3)*5; a,lam,r,m,t=rng.uniform(0.3,2),rng.choice([-2.,1.,3.]),rng.uniform(0.1,0.9),rng.uniform(0,2),rng.uniform(0,6)
    cs=rng.normal(size=3)+1j*rng.normal(size=3); chis=[rng.normal(size=4)+1j*rng.normal(size=4) for _ in ks]; th=rng.normal(size=2)+1j*rng.normal(size=2)
    lhs=0.; rhs=0.
    for ka,ca,ca_ in zip(ks,cs,chis):
        for kb,cb,cb_ in zip(ks,cs,chis):
            ph=np.conj(ca)*cb*np.exp(1j*(kb-ka)*t)
            lhs+=(ph*np.vdot(Psi(ca_,th),Psi(H(kb,a,lam,r,m)@cb_,th))).real
            Hbar=H((ka+kb)/2,a,lam,r,m); Xi=(ka+kb)/(2*a)*np.vdot(ca_,g5@cb_)+m*np.vdot(ca_,g0@Sig1@cb_)
            rhs+=(ph*(np.vdot(ca_,Hbar@cb_)*np.vdot(th,th)+Xi*np.vdot(th,s1@th))).real
    worst=max(worst,abs(lhs-rhs)/abs(lhs))
print('max relative |LHS-RHS| over 20 random draws (3 modes, lambda in {-2,1,3}, m in [0,2]):',f'{worst:.2e}')
