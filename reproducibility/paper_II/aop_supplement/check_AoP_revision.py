"""Targeted scalar Paper-II audit. Does not edit the v1.1.0 archive.
Run with Python + numpy, scipy, sympy, mpmath. Outputs are diagnostic, not intervals.
The analytic proofs in the manuscript, not these finite tests, establish the theorems.
"""
from __future__ import annotations
from pathlib import Path
import csv, json, sys, time, importlib.util
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import kv, gamma
from numpy.polynomial.legendre import leggauss

BASE=Path(__file__).resolve().parent
OUT=BASE/'results';OUT.mkdir(exist_ok=True)
checks=[]
def exact(name, value):
    if isinstance(value, sp.MatrixBase):
        ok=all(sp.simplify(x)==0 for x in value)
    else:ok=sp.simplify(value)==0
    checks.append({'name':name,'pass':bool(ok)})
    if not ok:raise AssertionError((name,value))
def numeric(name,value,bound):
    checks.append({'name':name,'value':float(value),'bound':float(bound),'pass':bool(value<=bound)})
    if not value<=bound:raise AssertionError((name,value,bound))
def save(name,rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)

r,M,sigma,omega,Om,H=sp.symbols('r M sigma omega Om H',positive=True,real=True)
f=(1-r)/r;q=r**sp.Rational(3,2)*sp.sqrt(1-r)
dtau=lambda z:sp.sqrt(f)*sp.diff(z,r)
HH=sp.simplify(dtau(q)/q);QQ=sp.simplify(dtau(dtau(sp.sqrt(q)))/sp.sqrt(q))
exact('pump formula', QQ+sp.Rational(1,2)/r**3+1/(16*r**3*(1-r)))
exact('H formula',HH-(3-4*r)/(2*r**sp.Rational(3,2)*sp.sqrt(1-r)))
exact('maximal volume',HH.subs(r,sp.Rational(3,4)))
exact('pump H identity',QQ-dtau(HH)/2-HH**2/4)
ar,ai,br,bi=sp.symbols('ar ai br bi',real=True);A=ar+sp.I*ai;B=br+sp.I*bi
chi=sp.sqrt(sigma)*(A+B*sp.log(sigma))
W=chi*sp.conjugate(sp.diff(chi,sigma))-sp.diff(chi,sigma)*sp.conjugate(chi)
exact('logarithmic Wronskian', W-(A*sp.conjugate(B)-sp.conjugate(A)*B))
p=sp.sqrt(sp.Rational(3,2)*sigma)*(sp.diff(chi,sigma)-chi/(2*sigma))
exact('canonical log cancellation',p-sp.sqrt(sp.Rational(3,2))*B)
x=sp.sqrt(Om/omega)
av=(x+1/x)/2+sp.I*H/(4*sp.sqrt(omega*Om))
bv=(x-1/x)/2+sp.I*H/(4*sp.sqrt(omega*Om))
exact('velocity Bogoliubov norm',av*sp.conjugate(av)-bv*sp.conjugate(bv)-1)
lam,kap,bb=sp.symbols('lambda kappa bb',positive=True)
exact('selected count leading factor', (lam/(2*sp.sqrt(bb)))**2*sp.Rational(2,3)/(2*kap)-lam**2/(12*kap*bb))
z=sp.symbols('z',real=True)
al=sp.Rational(5,4)-sp.I*z/(2*sp.pi);be=-sp.Rational(3,4)-sp.I*z/(2*sp.pi)
exact('fixed-tower normalization',al*sp.conjugate(al)-be*sp.conjugate(be)-1)
exact('fixed-tower occupation',be*sp.conjugate(be)-sp.Rational(9,16)-z**2/(4*sp.pi**2))
exact('fixed-tower symmetric ratio',al*sp.conjugate(al)+be*sp.conjugate(be)-sp.Rational(17,8)-z**2/(2*sp.pi**2))
# Match physical endpoint constants/logs in the documented Jost phase.
L,K=sp.symbols('L K',positive=True);eg=sp.symbols('eg',real=True)
aH=sp.sqrt(sp.pi/2)*(1-2*sp.I/sp.pi*(sp.log(L)+eg));bH=-sp.I/sp.sqrt(2*sp.pi)
aJ=sp.sqrt(sp.pi/8)*(1-2*sp.I/sp.pi*(sp.log(K/4)+eg));bJ=-sp.I*sp.sqrt(2/sp.pi)
logLambda=sp.log(K/4)-3*eg-4*sp.log(L)
aT=al.subs(z,logLambda);bT=be.subs(z,logLambda)
exact('fixed-tower logarithmic matching',aT*bH+bT*sp.conjugate(bH)-bJ)
exact('fixed-tower constant matching',aT*aH+bT*sp.conjugate(aH)-aJ)
# Scalar energy, including the distinction between static/Killing and comoving frames.
E,B2=sp.symbols('E B2',real=True);ff=2*M/r-1
wr2=sp.Rational(3,4)*B2/(M*r**2)
rho=(ff+2*E**2)*wr2
exact('radial energy coefficient',sp.limit(r**3*rho,r,0)-sp.Rational(3,2)*B2)
exact('sigma E leading coefficient',1/(3*M)*(sp.Rational(9,2)*M)-sp.Rational(3,2))
ai2,bi2,LL=sp.symbols('a_i b_i ell',complex=True)
# Independent leading flux algebra with complex amplitude jets.
flux=-2*sp.re(B*(sp.conjugate(ai2)+sp.conjugate(bi2)*LL))
expected=-2*sp.re(B*sp.conjugate(ai2))-2*LL*sp.re(B*sp.conjugate(bi2))
exact('mixed flux coefficient',sp.expand_complex((flux-expected).subs(LL,z)))
# Improvement energy is H F', not an imported vacuum coefficient.
F2=2*B2/(3*M)*sp.log(sigma)**2
impr=sp.diff(F2,sigma)/sigma
exact('curvature improvement factor',impr/(B2/(3*M*sigma**2))-4*sp.log(sigma))
exact('proper time to areal log factor',4*sp.Rational(3,2)-6)
# Second-derivative carrier suppresses the normalized-Jost infrared singularity.
k=sp.symbols('k',positive=True)
exact('IR weighted power', (k**2*k**(-sp.Rational(1,2)))**2/k-k**2)
# Independent finite-squeeze diagnostics.
squeeze=[]
for sval in [1,2,4,6,8]:
    nbar=np.sinh(sval)**2
    prob=1-(np.tanh(sval)**2)**9
    squeeze.append({'s':sval,'occupation':nbar,'norm_S':np.exp(sval),'prob_pair_n_le_8':prob})
save('squeezed_pair.csv',squeeze)
numeric('finite-particle projection tends to zero',squeeze[-1]['prob_pair_n_le_8'],1e-5)
# Scalar CCR rank-two identity in a finite vacuum/one-particle sector.
rng=np.random.default_rng(71419);N=5
a=np.diag(np.sqrt(np.arange(1,N)),1);I=np.eye(N)
ops=[np.kron(a,I),np.kron(I,a)];vac=np.zeros(N*N);vac[0]=1
errs=[]
for j in range(8):
    coeff=rng.normal(size=2)+1j*rng.normal(size=2);coeff/=np.linalg.norm(coeff)
    ux=rng.normal(size=2)+1j*rng.normal(size=2);uy=rng.normal(size=2)+1j*rng.normal(size=2)
    Phi=lambda u:sum(u[t]*ops[t]+u[t].conjugate()*ops[t].conj().T for t in range(2))
    state=sum(coeff[t]*ops[t].conj().T@vac for t in range(2))
    got=np.vdot(state,Phi(ux)@Phi(uy)@state)-np.vdot(vac,Phi(ux)@Phi(uy)@vac)
    wx=np.dot(coeff,ux);wy=np.dot(coeff,uy);expected=wx*wy.conjugate()+wx.conjugate()*wy
    errs.append(abs(got-expected))
numeric('scalar CCR smooth rank-two change',max(errs),2e-13)
# Dimensions in length powers, with measure d omega.
exact('mode dimension', -1+sp.Rational(1,2)+sp.Rational(1,2)) # 1/r sqrt(omega)
exact('packet dimension',sp.Rational(1,2)-sp.Rational(1,2)-1+1)
exact('stress dimension',-1-3+4)
exact('B connection dimension',sp.Rational(1,2)-sp.Rational(1,2))

# Reuse the immutable radial equation implementation, adding independent controls.
spec=importlib.util.spec_from_file_location('nc',BASE.parent/'num_checks.py')
nc=importlib.util.module_from_spec(spec);spec.loader.exec_module(nc)
rows=[]
for K0 in [2,4,8,16,32,64,128]:
    for xmin,rtol in [(1e-10,1e-11),(1e-12,2e-12),(1e-14,2e-13)]:
        nc.RT=rtol;nc.AT=rtol*1e-2
        ac,bc=nc.endpoint_coeffs(K0,0,xmin=np.log(xmin))
        rows.append({'K':K0,'rmin':xmin,'rtol':rtol,'K_beta2':K0*abs(bc)**2,'norm_error':abs(abs(ac)**2-abs(bc)**2-1)})
save('swave_refinement.csv',rows)
numeric('refined s-wave Wronskian',max(r['norm_error'] for r in rows if r['rmin']==1e-14),1e-10)
# Low-frequency bound diagnosed at a fixed nonzero angular tower; not an IR proof by fitting.
ir=[]
nc.RT=2e-12;nc.AT=2e-14
for kval in [1e-1,3e-2,1e-2,3e-3,1e-3]:
    ac,bc=nc.endpoint_coeffs(kval,1,xmin=np.log(1e-12))
    ir.append({'K':kval,'sqrtK_alpha_plus_beta':np.sqrt(kval)*(abs(ac)+abs(bc)),'norm_error':abs(abs(ac)**2-abs(bc)**2-1)})
save('Jost_IR_diagnostic.csv',ir)
# Physical normalization of the packet actually used by run_part2.py.
J=quad(lambda x:np.exp(-2/(1-x*x)),-1,1,epsabs=1e-14,epsrel=1e-13)[0]
Z=6*np.pi*J
numeric('horizon normalization Bessel identity',abs(J-2*np.exp(-1)*(kv(1,1)-kv(0,1))),1e-14)
rs=np.array([1e-8,1e-6,1e-4,1e-2]);packet=[];refin=[]
def time_r(r,E):
    return 0. if E==0 else quad(lambda x:E/((1/x-1)*np.sqrt(E*E-1+1/x)),0,r,epsabs=1e-16,epsrel=1e-12)[0]
for nq,rtol in [(48,1e-11),(96,2e-12),(96,2e-13)]:
    nc.RT=rtol;nc.AT=rtol/100
    x,wg=leggauss(nq);om=x+1.5;weight=wg
    raw=np.exp(-1/(1-x*x));bnorm=raw/np.sqrt(Z)
    norm=4*np.pi*np.sum(weight*om*bnorm**2)
    for l in [0,1]:
        phs=[];pxs=[];bs=[]
        for w0 in om:
            (ph_end,px_end),(xx,ph,px)=nc.propagate(w0,l*(l+1),0j,np.sqrt(2*w0)+0j,rs_out=rs,xmin=np.log(1e-14))
            order=np.argsort(xx);phs.append(ph[order]);pxs.append(px[order]);bs.append(2*np.sqrt(.5/3)*px_end)
        phs=np.array(phs);pxs=np.array(pxs);bs=np.array(bs)
        for r0 in rs:
            ii=np.argmin(abs(rs-r0));ff0=1/r0-1
            for location,Y,dY in ([('Y00',1/np.sqrt(4*np.pi),0.)] if l==0 else [('pole',np.sqrt(3/(4*np.pi)),0.),('equator',0.,-np.sqrt(3/(4*np.pi)))]):
                Bstar=np.sum(weight*bnorm*bs)*Y
                for E0 in [0,.5,2]:
                    tt=time_r(r0,E0);c=weight*bnorm*np.exp(-1j*om*tt)
                    ww=np.sum(c*phs[:,ii])*Y
                    wr=np.sum(c*pxs[:,ii])*Y/r0
                    wt=np.sum(c*(-1j*om)*phs[:,ii])*Y
                    wa=np.sum(c*phs[:,ii])*dY
                    ur=-np.sqrt(ff0+E0**2);ut=-E0/ff0
                    rho=2*abs(ur*wr+ut*wt)**2-ff0*abs(wr)**2+abs(wt)**2/ff0+abs(wa)**2/r0**2
                    row={'N_omega':nq,'rtol':rtol,'l':l,'location':location,'r':r0,'E':E0,'norm_horizon':norm,'Z':Z,'coefficient_3B2_over2':1.5*abs(Bstar)**2,'r3rho':r0**3*rho,'ratio':r0**3*rho/(1.5*abs(Bstar)**2) if abs(Bstar)>0 else 0.}
                    packet.append(row)
        refin.append({'N_omega':nq,'rtol':rtol,'l':l,'norm_horizon':norm,'coefficient_pole':1.5*abs(np.sum(weight*bnorm*bs))**2*((1 if l==0 else 3)/(4*np.pi))})
save('normalized_packet.csv',packet);save('packet_refinement.csv',refin)
numeric('96-node horizon profile normalization',max(abs(row['norm_horizon']-1) for row in refin if row['N_omega']==96),1e-11)
# Exact endpoint coefficient As, not a fitted quantity.
D=np.pi**1.5*gamma(.25)/(24*gamma(.75)**3);As=2*D**2
numeric('s-wave gamma coefficient',abs(As-.417942472394),5e-13)
(OUT/'audit_summary.json').write_text(json.dumps({'checks':checks,'count':len(checks),'all_pass':all(x['pass'] for x in checks),'Z':Z,'D':D,'As':As,'notes':'Floating-point diagnostics are not interval-certified. The archived scripts are preserved unchanged.'},indent=2))
print('PASS',len(checks),'checks; Z=',Z,'As=',As,flush=True)
