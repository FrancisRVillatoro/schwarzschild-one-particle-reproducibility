#!/usr/bin/env python3
"""Work4: independent CAR restriction, continuum normalization and Ward checks.

The continuum solution is a prescribed smooth spectral packet, NOT numerically
constructed global Unruh scattering data. The finite CAR model tests algebraic
restriction, not an astrophysical transmission coefficient.
"""
from pathlib import Path
import csv, json, platform, sys, time
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag
from numpy.polynomial.legendre import leggauss
from work3_stress_core import (G1,G2,g0,g5,S1,g,AA,harmonic_data,
                               tensor_modal,tensor_full,car_check,local_checks)
OUT=Path(__file__).resolve().parent
I4=np.eye(4,dtype=complex)


def covariance_restriction():
    """Explicit four-mode Fock representation with a non-projective restriction."""
    low=np.array([[0.,1.],[0.,0.]],complex); Z=np.diag([1.,-1.]); I2=np.eye(2)
    cs=[]
    for a in range(4):
        c=np.ones((1,1),complex)
        for b in range(4): c=np.kron(c,Z if b<a else low if b==a else I2)
        cs.append(c)
    vac=np.zeros(16,complex); vac[0]=1
    sea=cs[2].conj().T@cs[3].conj().T@vac
    P=np.diag([1.,1.,0.,0.])
    W=np.array([.6,.8j,0,0]); V=np.array([0,0,1.,-1j])/np.sqrt(2)
    # c^dagger(W) adds a particle; c(V) removes a filled negative mode.
    part=sum(W[a]*cs[a].conj().T@sea for a in range(4))
    anti=sum(V[a].conjugate()*cs[a]@sea for a in range(4))
    eta=np.array([.25,.70])
    J=np.vstack([np.diag(np.sqrt(eta)),np.diag(np.sqrt(1-eta))])
    ds=[sum(J[a,i].conjugate()*cs[a] for a in range(4)) for i in range(2)]
    def cov(st):
        return np.array([[np.vdot(st,ds[i]@ds[j].conj().T@st)
                          for j in range(2)] for i in range(2)])
    C=J.conj().T@P@J; w=J.conj().T@W; v=J.conj().T@V
    Cpart=C-np.outer(w,w.conj()); Canti=C+np.outer(v,v.conj())
    err=max(np.linalg.norm(cov(sea)-C), np.linalg.norm(cov(part)-Cpart),
            np.linalg.norm(cov(anti)-Canti))
    assert err<2e-14
    assert min(np.linalg.eigvalsh(Cpart))>-2e-14
    assert max(np.linalg.eigvalsh(Canti))<1+2e-14
    n=float(np.vdot(w,w).real)
    wrong=C-np.outer(w,w.conj())/n
    return dict(description='Finite CAR restriction; not a black-hole scattering model.',
                rank_one_error=float(err),interior_norm=n,
                restricted_covariance_eigenvalues=np.linalg.eigvalsh(C).tolist(),
                particle_covariance_eigenvalues=np.linalg.eigvalsh(Cpart).tolist(),
                antiparticle_covariance_eigenvalues=np.linalg.eigvalsh(Canti).tolist(),
                unjustified_renormalization_min_eigenvalue=float(min(np.linalg.eigvalsh(wrong))))


class Packet:
    """Massive/massless continuous diagnostic in regular x=sqrt(R/R_s)."""
    def __init__(self,N=96,Rs=2.3,mu=0.,l=2,rtol=2e-12):
        self.N=N;self.Rs=Rs;self.mu=mu;self.l=l
        x,w=leggauss(N);self.ks=-2+x;self.weights=w/(2*np.pi)
        amp=np.exp(-1/(1-x*x));amp/=np.sqrt(np.dot(self.weights,amp**2))
        # A genuine superposition of both transverse subspaces, with both
        # chiralities; normalization is deliberately N_int=0.37, not one.
        v=np.array([1.,.2+.4j,-.6+.3j,.7j]);v/=np.linalg.norm(v)
        self.b=np.sqrt(.37)*amp[:,None]*v
        self.Nint=float(np.dot(self.weights,np.sum(abs(self.b)**2,axis=1)))
        self.Kcharge=float(-np.dot(self.weights,self.ks*np.sum(abs(self.b)**2,axis=1)))
        self.A0=float(np.dot(self.weights,np.linalg.norm(self.b,axis=1)))
        self.A1=float(np.dot(self.weights,abs(self.ks)*np.linalg.norm(self.b,axis=1)))
        def rhs(x,Y):
            z=Y.reshape(N,4)
            h1=2*Rs*self.ks*x**3/(1-x*x)
            h2=2*l/np.sqrt(1-x*x)
            hm=2*Rs*mu*x*x/np.sqrt(1-x*x)
            return (1j*(h1[:,None]*(z@G1.T)+h2*(z@G2.T)+hm*(z@g0.T))).ravel()
        self.sol=solve_ivp(rhs,(0.,np.sqrt(.4)),self.b.ravel(),method='DOP853',
                           rtol=rtol,atol=rtol*1e-2,dense_output=True)
        if not self.sol.success: raise RuntimeError(self.sol.message)
    def z(self,R): return self.sol.sol(np.sqrt(R/self.Rs)).reshape(self.N,4)
    def tensor(self,X):
        R,T,theta,phi=X
        return tensor_modal(R,self.Rs,self.l,self.mu,theta,phi,self.ks,
                            self.z(R),self.weights,T)[0]
    def current(self,X):
        R,T,theta,phi=X; a=np.sqrt(self.Rs/R-1);q=a*R*R
        M,_,_,_,_=harmonic_data(self.l,theta,phi)
        psi=M@((self.weights*np.exp(1j*self.ks*T))@self.z(R))
        jf=np.array([np.vdot(psi,m@psi).real/q for m in
                     [I4,G1,g0@g[2],g0@g[3]]])
        return np.array([-a,1/a,1/R,1/(R*np.sin(theta))])*jf
    def mixed(self,X,drop_shear=False):
        R,T,theta,phi=X;a=np.sqrt(self.Rs/R-1)
        tf=self.tensor(X)
        if drop_shear:
            tf[1,2]=tf[2,1]=tf[1,3]=tf[3,1]=0.
        e=np.array([-1/a,a,R,R*np.sin(theta)])
        gi=np.array([a*a,-1/(a*a),-1/R**2,-1/(R*np.sin(theta))**2])
        return gi[:,None]*(e[:,None]*e[None,:]*tf)


def christoffel(R,Rs,theta):
    a2=Rs/R-1; ap=-Rs/R**2
    gd=np.array([1/a2,-a2,-R**2,-R**2*np.sin(theta)**2])
    dg=np.zeros((4,4))
    dg[0]=[-ap/a2**2,-ap,-2*R,-2*R*np.sin(theta)**2]
    dg[2,3]=-2*R**2*np.sin(theta)*np.cos(theta)
    C=np.zeros((4,4,4))
    for u in range(4):
        for v in range(4):
            for w in range(4):
                C[u,v,w]=((dg[v,u] if u==w else 0.)+(dg[w,u] if u==v else 0.)
                          -(dg[u,v] if v==w else 0.))/(2*gd[u])
    return C


def derivative(fun,X,j,h):
    e=np.zeros(4);e[j]=h
    return (fun(X-2*e)-8*fun(X-e)+8*fun(X+e)-fun(X+2*e))/(12*h)


def ward_test(packet):
    X=np.array([packet.Rs*.08,.31,1.07,.23])
    C=christoffel(X[0],packet.Rs,X[2]); tm=packet.mixed(X);jc=packet.current(X)
    trace=packet.tensor(X);R=X[0];a=np.sqrt(packet.Rs/R-1);q=a*R*R
    M,_,_,_,_=harmonic_data(packet.l,X[2],X[3]);psi=M@((packet.weights*np.exp(1j*packet.ks*X[1]))@packet.z(R))
    trace_error=abs(trace[0,0]-np.trace(trace[1:,1:])-packet.mu*np.vdot(psi,g0@psi).real/q)
    rows=[]
    for h in [8e-3,4e-3,2e-3,1e-3]:
        hs=h*np.array([X[0],packet.Rs,1.,1.])
        DT=np.array([derivative(packet.mixed,X,j,hs[j])[j,:] for j in range(4)])
        conn1=np.einsum('mml,ln->n',C,tm);conn2=np.einsum('lmn,ml->n',C,tm)
        res=DT.sum(axis=0)+conn1-conn2
        scale=np.sum(abs(DT),axis=0)+abs(conn1)+abs(conn2)
        terr=float(np.max(abs(res)/np.maximum(scale,1e-30)))
        DJ=np.array([derivative(packet.current,X,j,hs[j])[j] for j in range(4)])
        cJ=np.einsum('mml,l->',C,jc)
        jerr=float(abs(DJ.sum()+cJ)/(np.sum(abs(DJ))+abs(cJ)))
        rows.append(dict(mu=packet.mu,h=float(h),stress_relative_residual=terr,
                         current_relative_residual=jerr))
    assert rows[-1]['stress_relative_residual']<2e-8
    assert rows[-1]['current_relative_residual']<2e-8
    assert trace_error/max(1,np.max(abs(trace)))<5e-12
    return rows,float(trace_error)


def integrated_charges(p):
    # Direct quadratures in T and on the sphere, after evaluating the exact
    # angular Gram factors. They are independent of the spectral charge formula.
    xu,wu=leggauss(32)
    ints=np.array([harmonic_data(p.l,np.arccos(u),.27)[3:] for u in xu])
    s,d=(2*np.pi*wu)@ints
    rows=[]
    for Tmax,NT in [(30.,320),(60.,640),(90.,960)]:
        xt,wt=leggauss(NT);Ts=Tmax*xt;wt=Tmax*wt
        phase=np.exp(1j*Ts[:,None]*p.ks[None,:])*p.weights
        for rr in [.01,.1,.3]:
            R=p.Rs*rr;a=np.sqrt(p.Rs/R-1);q=a*R*R;z=p.z(R)
            Z=phase@z;Zk=phase@(p.ks[:,None]*z)
            quad=lambda A,Y:np.einsum('bi,ij,bj->b',Z.conj(),A,Y).real
            density=quad(s*I4+d*S1,Z)
            flux=-(quad(s*I4+d*S1,Zk)/a+p.l*d*quad(g5@G2,Z)/(2*R))/q
            num=float(wt@density);energy=float(q*a*(wt@flux))
            rows.append(dict(mu=p.mu,Tmax=Tmax,N_T=NT,r=rr,
                             number=num,number_exact=p.Nint,
                             number_error=abs(num-p.Nint),
                             Killing_charge=energy,Killing_charge_exact=p.Kcharge,
                             charge_error=abs(energy-p.Kcharge)))
    assert max(x['number_error'] for x in rows if x['Tmax']==90)<1e-8
    assert max(x['charge_error'] for x in rows if x['Tmax']==90)<1e-8
    return rows


def terminal_and_bound(p):
    T=.31;theta=1.07;phi=.23
    M,_,_,s,d=harmonic_data(p.l,theta,phi)
    B=(p.weights*np.exp(1j*p.ks*T))@p.b;psi0=M@B
    h=float(np.vdot(B,G2@B).real);limit=p.l*s*h/np.sqrt(p.Rs)
    shear=-3/8*np.vdot(psi0,1j*g[1]@g[2]@psi0).real
    rows=[]
    for r in [.01,.003,.001,3e-4,1e-4,1e-5,1e-6]:
        R=p.Rs*r;a=np.sqrt(1/r-1);q=a*R*R
        tensor=p.tensor(np.array([R,T,theta,phi]))
        Ipar=-R-p.Rs*np.log1p(-r)
        Sigma=p.Rs*(np.arcsin(np.sqrt(r))-np.sqrt(r*(1-r)))
        bound=((s+abs(d))/(q*a)*p.A0*p.A1
               +p.mu*(s+abs(d))/q*p.A0**2
               +2*p.l*s/(q*R)*p.A0*(Ipar*p.A1+p.mu*Sigma*p.A0))
        leading=p.l*s*h/(q*R)
        diff=float(abs(tensor[0,0]-leading))
        rows.append(dict(mu=p.mu,r=r,scaled_density=float(R**2.5*np.sqrt(1-r)*tensor[0,0]),
                         density_limit=limit,scaled_T12=float(R**3*tensor[1,2]),T12_limit=float(shear),
                         additive_error=diff,explicit_error_bound=float(bound),bound_ratio=diff/bound))
        assert diff<bound*(1+1e-7)
    return rows


def dimensions_and_box():
    # A smooth continuum profile has c_n=f(k_n)/sqrt(L_T), not f(k_n).
    def raw(k):
        x=k+2
        out=np.zeros_like(k);keep=abs(x)<1
        out[keep]=np.exp(-1/(1-x[keep]**2))
        return out
    x,w=leggauss(256);norm=np.sqrt(np.dot(w,raw(x-2)**2)/(2*np.pi))
    rows=[]
    for LT in [30.,60.,120.,240.]:
        dk=2*np.pi/LT;ks=dk*np.arange(np.ceil(-3/dk),np.floor(-1/dk)+1)
        c=raw(ks)/(norm*np.sqrt(LT))
        rows.append(dict(L_T=LT,delta_kappa=dk,box_norm=float(np.vdot(c,c).real),
                         number_of_modes=len(ks)))
    assert abs(rows[-1]['box_norm']-1)<2e-6
    # Full tensors, all 16 entries, under physical scale change.
    rng=np.random.default_rng(405);z=rng.normal(size=(3,4))+1j*rng.normal(size=(3,4))
    kh=np.array([-3.,-1.4,.2]);wh=np.array([.2,.3,.1]);ts=[]
    for Rs in [1.,2.3,7.]:
        t,_,_=tensor_modal(.13*Rs,Rs,2,.27/Rs,1.04,.2,kh/Rs,np.sqrt(Rs)*z,wh/Rs,.4*Rs)
        ts.append(Rs**4*t)
    err=max(np.linalg.norm(t-ts[0])/np.linalg.norm(ts[0]) for t in ts[1:])
    assert err<2e-13
    return rows,float(err)


def dumpcsv(name,rows):
    with (OUT/name).open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)


def main():
    start=time.time(); summary={}
    summary['CAR_restriction']=covariance_restriction()
    summary['CAR_particle_antiparticle']=car_check()
    summary['local_tensor_random']=local_checks()
    wr=[];cr=[];br=[];trace=[]
    for mu in [0.,.35]:
        p=Packet(mu=mu)
        rr,te=ward_test(p);wr+=rr;trace.append(te)
        cr+=integrated_charges(p);br+=terminal_and_bound(p)
    box,dim=dimensions_and_box()
    dumpcsv('ward_identity_refinement.csv',wr)
    dumpcsv('continuum_charges.csv',cr)
    dumpcsv('terminal_bounds.csv',br)
    dumpcsv('box_limit.csv',box)
    summary.update(max_trace_absolute_error=max(trace),physical_scale_relative_error=dim,
                   ward_rows=wr,charge_max_error_Tmax90=max(r['charge_error'] for r in cr if r['Tmax']==90),
                   number_max_error_Tmax90=max(r['number_error'] for r in cr if r['Tmax']==90),
                   max_terminal_error_to_bound=max(r['bound_ratio'] for r in br),
                   final_box_norm=box[-1]['box_norm'],elapsed_seconds=time.time()-start)
    summary['environment']={'Python':sys.version,'NumPy':np.__version__,'SciPy':scipy.__version__,'platform':platform.platform()}
    (OUT/'work4_results.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
