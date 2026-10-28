#!/usr/bin/env python3
"""Independent addition-kernel and blown-up-sphere checks (no OCR or external data)."""
import csv,json
from pathlib import Path
import numpy as np
from scipy.special import eval_jacobi,jv
from spin_geometry import harmonics,parallel,exact,xyz,s1,s2,s3,eye
OUT=Path(__file__).resolve().parent

def exp_point(theta,phi,Z):
 p=xyz(theta,phi);e1=np.array([np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),-np.sin(theta)])
 e2=np.array([-np.sin(phi),np.cos(phi),0.]);q=np.linalg.norm(Z)
 x=np.cos(q)*p+np.sinc(q/np.pi)*(Z[0]*e1+Z[1]*e2)
 return float(np.arccos(np.clip(x[2],-1,1))),float(np.arctan2(x[1],x[0]))

def blown(l,Z,Zp):
 theta,phi=1.1,.4;t,p=exp_point(theta,phi,np.array(Z)/l);tp,pp=exp_point(theta,phi,np.array(Zp)/l)
 A=parallel(theta,phi,t,p)[0]; B=parallel(tp,pp,theta,phi)[0]
 return (2*np.pi/l)*A@exact(l,t,p,tp,pp)@B

def flat(Z,Zp):
 V=np.array(Z)-Zp;z=np.linalg.norm(V)
 return .5*(jv(0,z)*eye-1j*(jv(1,z)/z if z else .5)*(V[0]*s2+V[1]*s3))

def run():
 errors=[];diagerr=[];res=[]
 pairs=[(.9,.2,.7,.2),(1.1,.7,.9,.1),(1.8,-.4,.8,.2),(1.4,.3,1.4,.3)]
 for l in [1,2,5,12,24]:
  for t,p,tp,pp in pairs:
   h=harmonics(l,t,p);hp=harmonics(l,tp,pp);K=h.T@hp.conj()
   errors.append(float(np.linalg.norm(K-exact(l,t,p,tp,pp))))
  h=harmonics(l,1.1,.3);diagerr.append(float(np.linalg.norm(h.T@h.conj()-l/(4*np.pi)*eye)))
  # Differential eigen-equation, computed from the exact invariant kernel.
  t,p,tp,pp=1.1,.7,.9,.1;d=2e-5
  dt=(-exact(l,t+2*d,p,tp,pp)+8*exact(l,t+d,p,tp,pp)-8*exact(l,t-d,p,tp,pp)+exact(l,t-2*d,p,tp,pp))/(12*d)
  dp=(-exact(l,t,p+2*d,tp,pp)+8*exact(l,t,p+d,tp,pp)-8*exact(l,t,p-d,tp,pp)+exact(l,t,p-2*d,tp,pp))/(12*d)
  K=exact(l,t,p,tp,pp);SK=s2@(dt+.5/np.tan(t)*K)+s3@dp/np.sin(t)
  res.append(float(np.linalg.norm(SK+1j*l*K)))
 rows=[];Z=np.array([.65,-.8]);Zp=np.array([-.35,.55])
 for l in [8,16,32,64,128]:
  E=blown(l,Z,Zp)-flat(Z,Zp);d=2e-4
  grad=(blown(l,Z+[d,0],Zp)-flat(Z+[d,0],Zp)-blown(l,Z-[d,0],Zp)+flat(Z-[d,0],Zp))/(2*d)
  rows.append(dict(lambda_=l,error=float(np.linalg.norm(E)),scaled_error=float(l*l*np.linalg.norm(E)),dZ1_error=float(np.linalg.norm(grad))))
 with (OUT/'invariant_harmonics_results.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 assert max(errors)<2e-12 and max(diagerr)<2e-12
 assert max(res)<2e-7
 for a,b in zip(rows,rows[1:]): assert 3.8<a['error']/b['error']<4.2
 summary={'spectral_sum_max_error':max(errors),'coincidence_max_error':max(diagerr),'differential_residual_max_error':max(res),'point_pairs':len(pairs),'multipoles':[1,2,5,12,24],'normal_chart':rows}
 (OUT/'invariant_harmonics_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 print(json.dumps(summary,indent=2))
if __name__=='__main__':run()
