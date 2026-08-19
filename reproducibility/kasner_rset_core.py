
import numpy as np, math
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss

def G0(s,L,u,mu):
    return (L*L*u*u*np.exp(8*s/3)
            +L*L*(1-u*u)*np.exp(2*s/3)
            +mu*mu*np.exp(2*s))

def G_derivs(s,L,u,mu):
    a=L*L*u*u*np.exp(8*s/3)
    b=L*L*(1-u*u)*np.exp(2*s/3)
    c=mu*mu*np.exp(2*s)
    aa=8/3; bb=2/3; cc=2
    return [a*aa**n+b*bb**n+c*cc**n for n in range(5)]

def Oders_from_F(F):
    f,f1,f2,f3,f4=F[:5]
    O=np.sqrt(f)
    O1=f1/(2*O)
    O2=(2*f*f2-f1*f1)/(4*f**1.5)
    O3=(4*f*f*f3-6*f*f1*f2+3*f1**3)/(8*f**2.5)
    O4=(-4*(4*f1*f3+3*f2*f2)*f*f
        +8*f**3*f4+36*f*f1*f1*f2-15*f1**4)/(16*f**3.5)
    return O,O1,O2,O3,O4

def W4_mode_s(s,L,u,mu):
    O,O1,O2,O3,O4=Oders_from_F(G_derivs(s,L,u,mu))
    w2=(-2*O*O2+3*O1**2)/(8*O**3)
    w4=-(-8*O**3*O4+80*O**2*O1*O3+52*O**2*O2**2
          -396*O*O1**2*O2+297*O1**4)/(128*O**7)
    return O+w2+w4

def start_s(L,u,mu,target=120.0):
    def Om(s): return np.sqrt(G0(s,L,u,mu))
    if Om(0)>=target: return 0.0
    hi=1.0
    while Om(hi)<target:
        hi*=1.5
    return brentq(lambda s:Om(s)-target,0,hi)

def exact_moments(L,u,mu,target=120.,rtol=1e-10):
    S=start_s(L,u,mu,target)
    W=W4_mode_s(S,L,u,mu)
    h=1e-5
    Wp=(W4_mode_s(S+h,L,u,mu)-W4_mode_s(S-h,L,u,mu))/(2*h)
    psi=1/np.sqrt(2*W)
    psip=(-Wp/(2*W)-1j*W)*psi
    if S>0:
        def rhs(s,z):
            p=z[0]+1j*z[1]; q=z[2]+1j*z[3]
            qq=-G0(s,L,u,mu)*p
            return [q.real,q.imag,qq.real,qq.imag]
        sol=solve_ivp(rhs,(S,0.0),[psi.real,psi.imag,psip.real,psip.imag],
                      method="DOP853",rtol=rtol,atol=rtol*1e-2)
        z=sol.y[:,-1]
        psi=z[0]+1j*z[1]; psip=z[2]+1j*z[3]
    wr=(psi*np.conj(psip)-psip*np.conj(psi)).imag
    sc=1/np.sqrt(wr)
    return abs(psi*sc)**2,abs(psip*sc)**2

def falling(a,n):
    y=1.
    for j in range(n): y*=a-j
    return y

def omega_derivs_sigma1(L,u,mu):
    F=[0.]*5
    F[0]=L*L+mu*mu
    for n in range(1,5):
        F[n]=L*L*(u*u*falling(2/3,n)+(1-u*u)*falling(-4/3,n))
    return Oders_from_F(F)

def ad4_moments(L,u,mu):
    O,O1,O2,O3,O4=omega_derivs_sigma1(L,u,mu)
    S,S1,S2,H=.25,-.5,1.5,1.
    w2=(4*O**2*S-2*O*O2+3*O1**2)/(8*O**3)
    w2p=(4*O**3*S1-4*O**2*O1*S-2*O**2*O3
         +10*O*O1*O2-9*O1**3)/(8*O**4)
    w4=-(16*O**4*S**2+16*O**4*S2-48*O**3*S*O2
         -80*O**3*O1*S1-8*O**3*O4+152*O**2*S*O1**2
         +80*O**2*O1*O3+52*O**2*O2**2
         -396*O*O1**2*O2+297*O1**4)/(128*O**7)
    A0=1/(2*O); A2=-w2/(2*O**2); A4=w2*w2/(2*O**3)-w4/(2*O**2)
    A=A0+A2+A4
    C0=O1/(2*O)+H/2
    C2=(w2p*O-O1*w2)/(2*O**2)
    D=O/2 + (w2/2+A0*C0*C0) + (w4/2+A2*C0*C0+2*A0*C0*C2)
    return A,D

def stress(A,D,L,u,mu):
    return np.array([
        .5*(D+(L*L+mu*mu)*A),
        .5*D+.5*(L*L*(2*u*u-1)-mu*mu)*A,
        .5*D-.5*(L*L*u*u+mu*mu)*A
    ])

def ren_mode(L,u,mu):
    A,D=exact_moments(L,u,mu)
    A4,D4=ad4_moments(L,u,mu)
    return stress(A,D,L,u,mu)-stress(A4,D4,L,u,mu)

def integrate_mu(mu,Nu=8,Nl=7):
    xu,wu=leggauss(Nu); us=(xu+1)/2; wu=wu/2
    cuts=sorted(set([0.0,max(8*mu,0.20),1.0,4.0,16.0,64.0]))
    cc=[cuts[0]]
    for c in cuts[1:]:
        if c>cc[-1]+1e-14: cc.append(c)
    total=np.zeros(3)
    pieces=[]
    for lo,hi in zip(cc[:-1],cc[1:]):
        xl,wl=leggauss(Nl)
        ls=(lo+hi)/2+(hi-lo)*xl/2; wl=wl*(hi-lo)/2
        piece=np.zeros(3)
        for L,wL in zip(ls,wl):
            av=np.zeros(3)
            for u,wU in zip(us,wu):
                av += wU*ren_mode(float(L),float(u),mu)
            piece += wL*L**2*av/(2*np.pi**2)
        total += piece
        pieces.append((lo,hi,*piece))
    return total,pieces


def integrate_mu_adapt(mu,Nu=10,Nl=8):
    xu,wu=leggauss(Nu); us=(xu+1)/2; wu=wu/2
    cuts=sorted(set([0.0,8*mu,0.20,1.0,4.0,16.0,64.0]))
    cc=[cuts[0]]
    for c in cuts[1:]:
        if c>cc[-1]+1e-14: cc.append(c)
    total=np.zeros(3); pieces=[]
    for lo,hi in zip(cc[:-1],cc[1:]):
        xl,wl=leggauss(Nl)
        ls=(lo+hi)/2+(hi-lo)*xl/2; wl=wl*(hi-lo)/2
        piece=np.zeros(3)
        for L,wL in zip(ls,wl):
            av=np.zeros(3)
            for u,wU in zip(us,wu):
                av += wU*ren_mode(float(L),float(u),mu)
            piece += wL*L**2*av/(2*np.pi**2)
        total += piece; pieces.append((lo,hi,*piece))
    return total,pieces


def ren_mode_massless_state(L,u,mu_reg):
    # Exact state is the massless Kasner positive-frequency state.
    A,D=exact_moments(L,u,0.0)
    A4,D4=ad4_moments(L,u,mu_reg)
    return stress(A,D,L,u,0.0)-stress(A4,D4,L,u,mu_reg)

def integrate_regulator(mu_reg,Nu=10,Nl=8):
    xu,wu=leggauss(Nu); us=(xu+1)/2; wu=wu/2
    cuts=sorted(set([0.0,8*mu_reg,0.20,1.0,4.0,16.0,64.0]))
    cc=[cuts[0]]
    for c in cuts[1:]:
        if c>cc[-1]+1e-14: cc.append(c)
    total=np.zeros(3); pieces=[]
    for lo,hi in zip(cc[:-1],cc[1:]):
        xl,wl=leggauss(Nl)
        ls=(lo+hi)/2+(hi-lo)*xl/2; wl=wl*(hi-lo)/2
        piece=np.zeros(3)
        for L,wL in zip(ls,wl):
            av=np.zeros(3)
            for u,wU in zip(us,wu):
                av += wU*ren_mode_massless_state(float(L),float(u),mu_reg)
            piece += wL*L**2*av/(2*np.pi**2)
        total += piece; pieces.append((lo,hi,*piece))
    return total,pieces


def ren_delta_mode(L,u,mu_reg):
    A,_=exact_moments(L,u,0.0)
    A4,_=ad4_moments(L,u,mu_reg)
    # p_parallel-p_perp for exact massless state minus regulator counterterm
    exact=.5*L*L*(3*u*u-1)*A
    sub=.5*L*L*(3*u*u-1)*A4
    return exact-sub

def integrate_delta(mu_reg,Nu=14,Nl=10):
    xu,wu=leggauss(Nu); us=(xu+1)/2; wu=wu/2
    cuts=sorted(set([0.0,8*mu_reg,0.20,1.0,4.0,16.0,64.0]))
    cc=[cuts[0]]
    for c in cuts[1:]:
        if c>cc[-1]+1e-14: cc.append(c)
    total=0.; pieces=[]
    for lo,hi in zip(cc[:-1],cc[1:]):
        xl,wl=leggauss(Nl)
        ls=(lo+hi)/2+(hi-lo)*xl/2; wl=wl*(hi-lo)/2
        piece=0.
        for L,wL in zip(ls,wl):
            av=0.
            for u,wU in zip(us,wu):
                av += wU*ren_delta_mode(float(L),float(u),mu_reg)
            piece += wL*L**2*av/(2*np.pi**2)
        total += piece; pieces.append((lo,hi,piece))
    return total,pieces


def integrate_delta_interval(mu_reg,lo,hi,Nu=22,Nl=16):
    xu,wu=leggauss(Nu); us=(xu+1)/2; wu=wu/2
    xl,wl=leggauss(Nl)
    ls=(lo+hi)/2+(hi-lo)*xl/2; wl=wl*(hi-lo)/2
    piece=0.
    for L,wL in zip(ls,wl):
        av=0.
        for u,wU in zip(us,wu):
            av += wU*ren_delta_mode(float(L),float(u),mu_reg)
        piece += wL*L**2*av/(2*np.pi**2)
    return piece


def delta_massless_mode(L,u):
    A,_=exact_moments(L,u,0.0)
    return .5*L*L*(3*u*u-1)*A

def integrate_delta_massless_interval(lo,hi,Nu=22,Nl=16):
    xu,wu=leggauss(Nu); us=(xu+1)/2; wu=wu/2
    xl,wl=leggauss(Nl)
    ls=(lo+hi)/2+(hi-lo)*xl/2; wl=wl*(hi-lo)/2
    piece=0.
    for L,wL in zip(ls,wl):
        av=0.
        for u,wU in zip(us,wu):
            av += wU*delta_massless_mode(float(L),float(u))
        # angularly integrated 4th-order counterterm:
        # int_0^1 (3u^2-1) A_ad^(4) du = -1/(63 L^3)
        # Since delta mode = .5 L^2 (...)*A and measure is L^2/(2pi^2),
        # the renormalized radial integrand is L^4/(4pi^2)
        # times [angular exact + 1/(63 L^3)].
        exact_angular = 2*av/(L*L)  # int(3u^2-1) A du
        radial = L**4/(4*np.pi**2)*(exact_angular + 1/(63*L**3))
        piece += wL*radial
    return piece
