"""Theorem F check (units 2M = 1, conventions of the manuscript).  Packet of 4-spinor modes with common lambda,
Delta rho = (a r^2)^{-1} [ A(t) |Theta|^2 + B(t) Theta^+ s1 Theta ],
A = Re sum conj(c_k) c_k' e^{i(k'-k)t} chi_k^+ Hbar_{kk'} chi_k',   B = Re sum conj(c_k) c_k' e^{i(k'-k)t} Xi_{kk'}.
Leading law:  r^{5/2}(1-r)^{1/2} Delta rho / (|lambda| |Theta|^2) -> q_+(t) - q_-(t),  q_pm = ||P_pm^perp sum_k c_k e^{ikt} chi_{I,k}(0)||^2,
and the Theta^+ s1 Theta coefficient r^{5/2}(1-r)^{1/2} B/|lambda| -> 0 like r^{3/2}.  Modes started in the instantaneous positive
eigenspace at r0 = 0.9 with chirality gamma^5 = +1 (Weyl packet); the t-average of the leading coefficient is sum |c_k|^2 (1 - 2 n_perp(k))."""
import numpy as np
from dirac_phase0 import H, proj, g0, g2, evolve
I2=np.eye(2); Z=np.zeros((2,2)); s1=np.array([[0,1],[1,0]],complex)
g5=np.block([[Z,I2],[I2,Z]]).astype(complex); Sig1=np.block([[s1,Z],[Z,s1]])
lam=1.; m=0.; ks=[200,210]; cs=[1/np.sqrt(2),1/np.sqrt(2)]; rs=np.array([1e-2,1e-3,1e-4,1e-5,1e-6,1e-7])
# chirality-fixed modes (gamma^5 = +1): the positive eigenspace of H is two-dimensional and commutes with gamma^5 for m = 0
def evolve_chiral(k,lam,r0,rs):
    from scipy.integrate import solve_ivp
    w,v=np.linalg.eigh(H(r0,k,lam,0.)); e=v[:,w>0]; Pc=(np.eye(4)+g5)/2
    u0=Pc@e[:,0]
    if np.linalg.norm(u0)<0.5: u0=Pc@e[:,1]
    u0=u0/np.linalg.norm(u0)
    rhs=lambda r,y: np.concatenate([(-1j*(H(r,k,lam,0.)/np.sqrt(1/r-1))@(y[:4]+1j*y[4:])).real,(-1j*(H(r,k,lam,0.)/np.sqrt(1/r-1))@(y[:4]+1j*y[4:])).imag])
    sol=solve_ivp(rhs,(r0,rs[-1]),np.concatenate([u0.real,u0.imag]),t_eval=rs,method='DOP853',rtol=1e-11,atol=1e-13)
    return [sol.y[:4,i]+1j*sol.y[4:,i] for i in range(len(rs))]
data={k:evolve_chiral(k,lam,0.9,rs) for k in ks}
Pp=proj(g0@g2*np.sign(lam),+1); Pm=proj(g0@g2*np.sign(lam),-1)
chi0={k:data[k][-1] for k in ks}; nk={k:np.vdot(chi0[k],Pm@chi0[k]).real for k in ks}
print('terminal occupations n_perp:',{k:round(v,5) for k,v in nk.items()},'   sum |c|^2 (1-2n) =',round(sum(c**2*(1-2*nk[k]) for k,c in zip(ks,cs)),5))
ts=np.linspace(0,2*np.pi/10,13)[:-1]; leads=[]
print('     t      q_+-q_-     A r^{5/2}(1-r)^{1/2}/|lambda| at r=1e-2 .. 1e-7                               B r^{5/2}(1-r)^{1/2}/|lambda| at r=1e-2, 1e-3, 1e-4')
for t in ts:
    S0=sum(c*np.exp(1j*k*t)*chi0[k] for k,c in zip(ks,cs)); lead=np.vdot(S0,Pp@S0).real-np.vdot(S0,Pm@S0).real; leads.append(lead)
    rowA=[]; rowB=[]
    for i,r in enumerate(rs):
        a=np.sqrt(1/r-1); A=0.; B=0.
        for ka,ca in zip(ks,cs):
            for kb,cb in zip(ks,cs):
                ph=np.conj(ca)*cb*np.exp(1j*(kb-ka)*t)
                A+=(ph*np.vdot(data[ka][i],H(r,(ka+kb)/2,lam,m)@data[kb][i])).real
                B+=(ph*((ka+kb)/(2*a)*np.vdot(data[ka][i],g5@data[kb][i])+m*np.vdot(data[ka][i],g0@Sig1@data[kb][i]))).real
        rowA.append(A*r**2.5*np.sqrt(1-r)/lam/(a*r**2)); rowB.append(B*r**2.5*np.sqrt(1-r)/lam/(a*r**2))
    print(f'  {t:6.3f}   {lead:+.5f}   '+'  '.join(f'{v:+.5f}' for v in rowA)+'     '+'  '.join(f'{v:+.2e}' for v in rowB[:3]))
print('t-average of the leading coefficient over one beat period:',round(float(np.mean(leads)),5),'  min over t:',round(float(np.min(leads)),5))
