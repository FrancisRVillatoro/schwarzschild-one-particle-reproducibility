#!/usr/bin/env python3
"""CAR, tensor and continuous-packet checks for Work3.

The packet experiment uses explicitly prescribed continuous Cauchy/terminal
coefficients. It is NOT an exterior-scattering computation of an Unruh packet.
Only numpy and scipy are required.
"""
import csv,json
from pathlib import Path
import numpy as np
from scipy.linalg import block_diag,expm
from scipy.integrate import solve_ivp,simpson
from numpy.polynomial.legendre import leggauss
from spin_geometry import s1,s2,s3,harmonics,up,um
OUT=Path(__file__).resolve().parent
I2=np.eye(2);I4=np.eye(4);O2=np.zeros((2,2))
g0=block_diag(I2,-I2)
g=[g0]+[np.block([[O2,s],[-s,O2]]) for s in [s1,s2,s3]]
G1=g0@g[1];G2=np.block([[O2,-s3],[-s3,O2]]);g5=np.block([[O2,I2],[I2,O2]])
S1=block_diag(s1,s1);AA=[I4]+[-g0@x for x in g[1:]]

def harmonic_data(l,theta,phi,idx=0):
 j=l-.5;ms=np.arange(j,-j-1,-1);N=len(ms);jp=np.zeros((N,N),complex)
 for b in range(1,N):jp[b-1,b]=np.sqrt((j-ms[b])*(j+ms[b]+1))
 jy=(jp-jp.conj().T)/(2j);d=expm(-1j*theta*jy);dt=-1j*jy@d
 ca=np.flatnonzero(np.isclose(ms,.5))[0];cb=np.flatnonzero(np.isclose(ms,-.5))[0]
 pre=np.sqrt(2*l/(8*np.pi))*np.exp(1j*ms[idx]*phi)
 th=pre*(d[idx,ca]*up+d[idx,cb]*um);dh=pre*(dt[idx,ca]*up+dt[idx,cb]*um)
 ph=1j*ms[idx]*th
 T=np.column_stack([th,s1@th]);dT=np.column_stack([dh,s1@dh]);pT=np.column_stack([ph,s1@ph])
 return block_diag(T,T),block_diag(dT,dT),block_diag(pT,pT),float(np.vdot(th,th).real),float(np.vdot(th,s1@th).real)

def tensor_full(R,Rs,theta,Psi,dSigma,dT,dTheta,dPhi):
 a=np.sqrt(Rs/R-1);q=a*R*R;adot=-Rs/(2*R*R);rdot=a
 C=[np.zeros((4,4),complex),-adot/(2*a)*G1,-rdot/(2*R)*(g0@g[2]),-rdot/(2*R)*(g0@g[3])+1/(2*R*np.tan(theta))*(g[2]@g[3])]
 Ds=[-dSigma,dT/a,dTheta/R,dPhi/(R*np.sin(theta))]
 D=[Ds[i]+C[i]@Psi for i in range(4)]
 T=np.empty((4,4))
 for i in range(4):
  for j in range(4):
   # Each AA is Hermitian. Form the real antisymmetrized bilinear
   # without subtracting large nearly equal complex real parts.
   z1=np.vdot(Psi,AA[i]@D[j]); z2=np.vdot(Psi,AA[j]@D[i])
   z1b=np.vdot(D[j],AA[i]@Psi); z2b=np.vdot(D[i],AA[j]@Psi)
   scale=max(1.,abs(z1),abs(z2),abs(z1b),abs(z2b))
   assert max(abs(z1b-z1.conjugate()),abs(z2b-z2.conjugate()))<1e-12*scale
   T[i,j]=-(z1.imag+z2.imag)/(2*q)
 return T

def tensor_modal(R,Rs,l,mu,theta,phi,ks,z,weights=None,Tcoord=0.):
 M,Mt,Mp,s,d=harmonic_data(l,theta,phi)
 a=np.sqrt(Rs/R-1);q=a*R*R
 if weights is None:weights=np.ones(len(ks))
 e=weights*np.exp(1j*ks*Tcoord);F=e@z;Fk=(e*ks)@z
 Hz=np.array([((k/a)*G1+(l/R)*G2+mu*g0)@v for k,v in zip(ks,z)])
 HS=e@Hz
 tf=tensor_full(R,Rs,theta,M@F,1j*M@HS,1j*M@Fk,Mt@F,Mp@F)
 modal00=(s*(np.vdot(F,G1@Fk).real/a+l*np.vdot(F,G2@F).real/R+mu*np.vdot(F,g0@F).real)+d*(np.vdot(F,g5@Fk).real/a+mu*np.vdot(F,g0@S1@F).real))/q
 modal11=np.vdot(F,(s*G1+d*g5)@Fk).real/(q*a)
 modal01=-(np.vdot(F,(s*I4+d*S1)@Fk).real/a+l*d*np.vdot(F,g5@G2@F).real/(2*R))/q
 return tf,np.array([modal00,modal01,modal11]),M@F

def car_check():
 low=np.array([[0,1],[0,0]],complex);Z=np.diag([1,-1]);ops=[]
 for j in range(4):
  A=np.ones((1,1),complex)
  for k in range(4): A=np.kron(A,Z if k<j else low if k==j else I2)
  ops.append(A)
 vac=np.zeros(16,complex);vac[0]=1
 rng=np.random.default_rng(1021)
 U=rng.normal(size=(4,2))+1j*rng.normal(size=(4,2));V=rng.normal(size=(4,2))+1j*rng.normal(size=(4,2))
 f=np.array([.6,.8j]);h=np.array([1.,-1j])/np.sqrt(2)
 sp=sum(f[a]*ops[a].conj().T@vac for a in range(2));sm=sum(h[a]*ops[a+2].conj().T@vac for a in range(2))
 psi=[sum(U[i,a]*ops[a]+V[i,a]*ops[a+2].conj().T for a in range(2)) for i in range(4)]
 def cov(state): return np.array([[np.vdot(state,psi[i]@psi[j].conj().T@state) for j in range(4)] for i in range(4)])
 dp=cov(sp)-cov(vac);dm=cov(sm)-cov(vac)
 w=U@f;v=V@h.conj()  # b^dagger(h) yields v coefficients conjugated in <b b^dagger>.
 ep=np.linalg.norm(dp+np.outer(w,w.conj()));em=np.linalg.norm(dm-np.outer(v,v.conj()))
 assert ep<2e-13 and em<2e-13
 return dict(particle_rank_one_error=float(ep),antiparticle_rank_one_error=float(em),particle_state_norm=float(np.vdot(sp,sp).real),antiparticle_state_norm=float(np.vdot(sm,sm).real))

def local_checks():
 rng=np.random.default_rng(146);errs=[];shearerr=[];traceerr=[]
 for _ in range(12):
  R=float(rng.uniform(.08,.7));Rs=1.7;l=int(rng.integers(1,5));mu=float(rng.random());th=float(rng.uniform(.6,2.1));ph=float(rng.uniform(-.6,.6));ks=np.array([-3.2,-.4,1.7]);z=rng.normal(size=(3,4))+1j*rng.normal(size=(3,4))
  tf,tm,psi=tensor_modal(R,Rs,l,mu,th,ph,ks,z)
  exacts=np.array([tf[0,0],tf[0,1],tf[1,1]])
  errs.append(float(np.max(abs(exacts-tm))/max(1,np.max(abs(tm)))))
  q=np.sqrt(Rs/R-1)*R*R
  trace=tf[0,0]-tf[1,1]-tf[2,2]-tf[3,3]
  predicted=mu*np.vdot(psi,g0@psi).real/q
  traceerr.append(float(abs(trace-predicted)/max(1.,np.max(abs(tf)),abs(predicted))))
  # Isolate connection-only tensor, compare the analytical mixed-spatial coefficient.
  a=np.sqrt(Rs/R-1);q=a*R*R;bi=[None,(-Rs/(2*R*R))/a,a/R,a/R]
  C=[np.zeros((4,4),complex)]+[-bi[i]/2*g0@g[i] for i in range(1,4)]
  for A in [2,3]:
   mat=1j/(4*q)*(AA[1]@C[A]-C[A].conj().T@AA[1]+AA[A]@C[1]-C[1].conj().T@AA[A])
   pred=-(3*Rs-2*R)/(8*R**3*(Rs-R))*(1j*g[1]@g[A])
   shearerr.append(float(np.linalg.norm(mat-pred)))
 assert max(errs)<2e-13 and max(shearerr)<2e-11 and max(traceerr)<2e-13
 return dict(local_tensor_relative_error=max(errs),spin_shear_matrix_error=max(shearerr),random_packets=12,on_shell_trace_relative_error=max(traceerr))

def continuum_check(N=96,rtol=2e-12,atol=2e-14,save=True):
 Rs=3.;l=1;theta=1.05;phi=.3;Tcoord=.37
 x,w=leggauss(N);ks=-5+2*x;weights=w*2/(2*np.pi)
 amp=np.exp(-1/(1-x*x));amp/=np.sqrt(np.dot(weights,amp*amp))
 u=np.array([np.cos(1.),np.exp(.4j)*np.sin(1.)]);chi0=np.r_[u,u]/np.sqrt(2)
 z0=amp[:,None]*chi0
 def rhs(x,Y):
  z=Y.reshape(N,4)
  B=2*Rs*ks*x**3/(1-x*x);C=2*l/np.sqrt(1-x*x)
  return (1j*(B[:,None]*(z@G1.T)+C*(z@G2.T))).ravel()
 rhos=np.array([1e-2,3e-3,1e-3,3e-4,1e-4,3e-5,1e-5,1e-6])
 sol=solve_ivp(rhs,(0,np.sqrt(max(rhos))),z0.ravel(),method='DOP853',rtol=rtol,atol=atol,dense_output=True)
 if not sol.success:raise RuntimeError(sol.message)
 M,Mt,Mp,s,d=harmonic_data(l,theta,phi)
 B=(weights*np.exp(1j*ks*Tcoord))@z0;psi0=M@B
 h=np.vdot(B,G2@B).real;coef=l*s*h/np.sqrt(Rs)
 shearcoef=-3/8*np.vdot(psi0,1j*g[1]@g[2]@psi0).real
 rows=[];normerrs=[]
 for r in rhos:
  R=Rs*r;z=sol.sol(np.sqrt(r)).reshape(N,4);tf,tm,psi=tensor_modal(R,Rs,l,0,theta,phi,ks,z,weights,Tcoord)
  norm=np.dot(weights,np.sum(abs(z)**2,axis=1));normerrs.append(abs(norm-1))
  rows.append(dict(r=float(r),R=float(R),scaled_rho=float(R**2.5*np.sqrt(1-r)*tf[0,0]),rho_limit=float(coef),scaled_shear=float(R**3*tf[1,2]),shear_limit=float(shearcoef),interior_norm=float(norm)))
 if save:
  with (OUT/'continuum_packet_results.csv').open('w') as f:
   wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
 assert max(normerrs)<2e-10
 assert abs(rows[-1]['scaled_rho']/coef-1)<1e-6
 assert abs(rows[-1]['scaled_shear']/shearcoef-1)<.01
 # Dimension check at fixed dimensionless R, k, t and appropriately scaled f.
 r=.023;ks_hat=np.array([-5.,-2.]);z_hat=np.array([chi0,chi0*np.exp(.4j)]);wh=np.array([.11,.19]);th=.31
 tensors=[]
 for rs in [1.,3.,7.]:
  tf,_,_=tensor_modal(rs*r,rs,l,0,theta,phi,ks_hat/rs,np.sqrt(rs)*z_hat,wh/rs,rs*th)
  tensors.append(tf*rs**4)
 dimerr=max(float(np.linalg.norm(t-tensors[0])/np.linalg.norm(tensors[0])) for t in tensors[1:]);assert dimerr<1e-13
 return dict(description='Continuous terminal/Cauchy diagnostic, not an Unruh transmission calculation.',quadrature_points=N,max_norm_error=float(max(normerrs)),rho_leading_coefficient=float(coef),T12_leading_coefficient=float(shearcoef),mass_scale_relative_error=dimerr,rows=rows)

def run():
 out={'CAR':car_check(),'bilinear':local_checks(),'continuum':continuum_check()}
 (OUT/'continuum_stress_summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':run()
