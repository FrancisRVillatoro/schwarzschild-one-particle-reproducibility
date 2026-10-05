"""Convergence study for C_dyn and Delta_dyn (Appendix B of paper I).
Baseline parameters are those of reproduce_kasner_phi2.py and kasner_complete_rset.py."""
import numpy as np, sys, json, time
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
src=open('reproduce_kasner_phi2.py').read().split('raw=sum(')[0]; exec(src)   # Fd, W2d, amp2, g, segment
import kasner_rset_core as K

# ------------------------- C_dyn -------------------------
def Cdyn(n=16,Nu=24,target=60.,rtol=3e-9,cut=8.,fitL=(6.,8.,10.,12.)):
    segs=[(0,.25),(.25,1),(1,3),(3,cut)]
    raw=sum(segment(a,b,n,Nu,target,rtol) for a,b in segs)
    c7=-1/2835; Ls=np.array(fitL)
    gs=np.array([g(L,Nu=28,target=250.,rtol=3e-10) for L in Ls])
    res=gs-c7*Ls**-5; M=np.column_stack([Ls**-7,Ls**-9]); c9,c11=np.linalg.lstsq(M,res,rcond=None)[0]
    tail=(c7/(4*cut**4)+c9/(6*cut**6)+c11/(8*cut**8))
    return (raw+tail)/(2*np.pi**2), tail/(2*np.pi**2)

# ------------------------- Delta_dyn -------------------------
def exact_A(L,u,target,rtol):
    S=K.start_s(L,u,0.0,target); W=K.W4_mode_s(S,L,u,0.0); h=1e-5
    Wp=(K.W4_mode_s(S+h,L,u,0.0)-K.W4_mode_s(S-h,L,u,0.0))/(2*h)
    psi=1/np.sqrt(2*W); psip=(-Wp/(2*W)-1j*W)*psi
    if S>0:
        def rhs(s,z):
            p=z[0]+1j*z[1]; q=z[2]+1j*z[3]; qq=-K.G0(s,L,u,0.0)*p
            return [q.real,q.imag,qq.real,qq.imag]
        sol=solve_ivp(rhs,(S,0.0),[psi.real,psi.imag,psip.real,psip.imag],method="DOP853",rtol=rtol,atol=rtol*1e-2)
        z=sol.y[:,-1]; psi=z[0]+1j*z[1]; psip=z[2]+1j*z[3]
    wr=(psi*np.conj(psip)-psip*np.conj(psi)).imag
    return abs(psi)**2/wr
def Dpiece(lo,hi,Nu,Nl,target,rtol):
    xu,wu=leggauss(Nu); us=(xu+1)/2; wu=wu/2; xl,wl=leggauss(Nl)
    ls=(lo+hi)/2+(hi-lo)*xl/2; wl=wl*(hi-lo)/2; piece=0.
    for L,wL in zip(ls,wl):
        ang=sum(wU*(3*u*u-1)*exact_A(float(L),float(u),target,rtol) for u,wU in zip(us,wu))
        piece+=wL*L**4/(4*np.pi**2)*(ang+1/(63*L**3))
    return piece
def Ddyn(Nu=26,Nl=18,target=120.,rtol=1e-10,cut=16.,use_d9=True):
    edges=[0.0,0.2,1.0,4.0,cut] if cut>4 else [0.0,0.2,1.0,cut]
    pieces=[Dpiece(a,b,Nu,Nl,target,rtol) for a,b in zip(edges[:-1],edges[1:])]
    d7=4/2835; a,b=edges[-2],edges[-1]
    lead=d7/(8*np.pi**2)*(a**-2-b**-2); cd9=(a**-4-b**-4)/(16*np.pi**2)
    d9=(pieces[-1]-lead)/cd9 if use_d9 else 0.0
    tail=d7/(8*np.pi**2*b**2)+d9/(16*np.pi**2*b**4)
    return sum(pieces)+tail, tail, d9

which=sys.argv[1]; names=sys.argv[2:]
RC={'baseline':{}, 'target120':{'target':120.}, 'Nu32':{'Nu':32}, 'n24':{'n':24}, 'rtol3e-11':{'rtol':3e-11},
    'cut6':{'cut':6.}, 'fit8-14':{'fitL':(8.,10.,12.,14.)}, 'Nu40':{'Nu':40}, 'target240':{'target':240.}}
RD={'baseline':{}, 'target240':{'target':240.}, 'Nu34':{'Nu':34}, 'Nl26':{'Nl':26}, 'rtol1e-12':{'rtol':1e-12},
    'cut8':{'cut':8.}, 'no_d9':{'use_d9':False}}
t0=time.time()
for k in names:
    if which=='C': v=Cdyn(**RC[k])
    else: v=Ddyn(**RD[k])
    with open(f'conv_{which}.txt','a') as f: f.write(k+' '+' '.join(repr(float(x)) for x in v)+'\n')
    print(k,v,f'{time.time()-t0:.0f}s',flush=True)
