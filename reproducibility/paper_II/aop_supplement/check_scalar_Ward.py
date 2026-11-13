"""Independent scalar stress identity at random on-shell two-jets.
This tests the tensor algebra and coordinates, not a global packet or a theorem.
"""
from pathlib import Path
import json
import numpy as np
import sympy as sp
r,t,th,ph,M=sp.symbols('r t th ph M',positive=True,real=True)
coords=[r,t,th,ph];f=2*M/r-1
g=sp.diag(-1/f,f,r*r,r*r*sp.sin(th)**2);gi=g.inv()
dg=[g.diff(x) for x in coords];dgi=[gi.diff(x) for x in coords]
Gam=sp.MutableDenseNDimArray.zeros(4,4,4)
for a in range(4):
 for b in range(4):
  for c in range(4):
   Gam[a,b,c]=sp.simplify(sum(gi[a,d]*(sp.diff(g[d,c],coords[b])+sp.diff(g[d,b],coords[c])-sp.diff(g[b,c],coords[d]))/2 for d in range(4)))
gfun=sp.lambdify((r,th,M),g,'numpy');gifun=sp.lambdify((r,th,M),gi,'numpy')
dgfun=[sp.lambdify((r,th,M),x,'numpy') for x in dg]
dgifun=[sp.lambdify((r,th,M),x,'numpy') for x in dgi]
Gfun=sp.lambdify((r,th,M),[[[Gam[a,b,c] for c in range(4)] for b in range(4)] for a in range(4)],'numpy')
rng=np.random.default_rng(20261019);rows=[]
for case in range(20):
 m0=rng.uniform(.4,3.);r0=2*m0*rng.uniform(.08,.8);theta=rng.uniform(.2,2.8)
 G=np.array(Gfun(r0,theta,m0));gg=gfun(r0,theta,m0);iv=gifun(r0,theta,m0)
 dgg=np.array([a(r0,theta,m0) for a in dgfun]);div=np.array([a(r0,theta,m0) for a in dgifun])
 u=rng.normal(size=4)+1j*rng.normal(size=4)
 h=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));h=(h+h.T)/2
 rhs=sum(iv[b,c]*G[a,b,c]*u[a] for a in range(4) for b in range(4) for c in range(4))
 h[0,0]=(rhs-sum(iv[b,c]*h[b,c] for b in range(4) for c in range(4) if (b,c)!=(0,0)))/iv[0,0]
 gradnorm=np.einsum('ab,a,b',iv,u,u.conj()).real
 T=2*np.real(np.outer(u,u.conj()))-gg*gradnorm
 dT=[]
 for l in range(4):
  ds=(np.einsum('ab,a,b',div[l],u,u.conj())+2*np.real(np.einsum('ab,a,b',iv,h[l],u.conj()))).real
  dT.append(2*np.real(np.outer(h[l],u.conj())+np.outer(u,h[l].conj()))-dgg[l]*gradnorm-gg*ds)
 dT=np.array(dT);mix=iv@T;dmix=np.array([div[l]@T+iv@dT[l] for l in range(4)])
 residual=np.zeros(4)
 norm=np.zeros(4)
 for n in range(4):
  for mu in range(4):
   residual[n]+=dmix[mu,mu,n];norm[n]+=abs(dmix[mu,mu,n])
   for rho in range(4):
    a=G[mu,mu,rho]*mix[rho,n];b=G[rho,mu,n]*mix[mu,rho]
    residual[n]+=a-b;norm[n]+=abs(a)+abs(b)
 err=max(abs(residual)/np.maximum(norm,1))
 trerr=abs(np.trace(mix)+2*gradnorm)/(1+abs(gradnorm))
 rows.append({'case':case,'r':r0,'M':m0,'Ward_relative':float(err),'trace_relative':float(trerr)})
assert max(x['Ward_relative'] for x in rows)<1e-13
assert max(x['trace_relative'] for x in rows)<1e-13
out=Path(__file__).parent/'results'/'scalar_Ward_jets.json'
out.write_text(json.dumps({'all_pass':True,'kind':'random on-shell jet identities, not a global wavepacket','cases':rows},indent=2))
print('PASS 20 on-shell scalar stress jets; max divergence',max(x['Ward_relative'] for x in rows),'max trace',max(x['trace_relative'] for x in rows))
