import numpy as np
from scipy.linalg import expm
from scipy.special import eval_jacobi, jv
from scipy.integrate import quad
s1=np.array([[0,1],[1,0]],complex);s2=np.array([[0,-1j],[1j,0]],complex);s3=np.diag([1.,-1.]).astype(complex);eye=np.eye(2)
up=np.array([1.,1.])/np.sqrt(2);um=np.array([1.,-1.])/np.sqrt(2)
def harmonics(l,theta,phi):
 j=l-.5;m=np.arange(j,-j-1,-1);N=len(m);jp=np.zeros((N,N),complex)
 for b in range(1,N):
  v=m[b];jp[b-1,b]=np.sqrt((j-v)*(j+v+1))
 jy=(jp-jp.conj().T)/(2j)
 d=expm(-1j*theta*jy).real
 a=np.flatnonzero(np.isclose(m,.5))[0];b=np.flatnonzero(np.isclose(m,-.5))[0]
 return np.sqrt(2*l/(8*np.pi))*np.exp(1j*m[:,None]*phi)*(d[:,a,None]*up+d[:,b,None]*um)
def xyz(theta,phi):return np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])
def parallel(t,p,tp,pp):
 x=xyz(t,p);xp=xyz(tp,pp);g=np.arctan2(np.linalg.norm(np.cross(xp,x)),xp@x)
 if g<1e-12:return eye.copy(), np.array([1.,0.]),g
 v=(x-np.cos(g)*xp)/np.sin(g)
 def it(s):
  z=np.cos(s)*xp+np.sin(s)*v;dz=-np.sin(s)*xp+np.cos(s)*v
  return z[2]*(z[0]*dz[1]-z[1]*dz[0])/(z[0]**2+z[1]**2)
 beta=quad(it,0,g,epsabs=1e-13,epsrel=1e-13)[0]
 U=expm(.5j*beta*s1)
 tang=(np.cos(g)*x-xp)/np.sin(g)
 et=np.array([np.cos(t)*np.cos(p),np.cos(t)*np.sin(p),-np.sin(t)])
 ep=np.array([-np.sin(p),np.cos(p),0.])
 return U,np.array([tang@et,tang@ep]),g
def exact(l,t,p,tp,pp):
 U,n,g=parallel(t,p,tp,pp)
 f=np.cos(g/2)*eval_jacobi(l-1,0,1,np.cos(g));gg=np.sin(g/2)*eval_jacobi(l-1,1,0,np.cos(g))
 return l/(4*np.pi)*(f*eye-1j*gg*(n[0]*s2+n[1]*s3))@U
if __name__=='__main__':
 for l in [1,2,5,12]:
  for t,p,tp,pp in [(.9,.2,.7,.2),(1.1,.7,.9,.1),(1.8,-.4,.8,.2)]:
   h=harmonics(l,t,p);hp=harmonics(l,tp,pp)
   K=h.T@hp.conj()
   print(l,t,p,np.max(abs(K-exact(l,t,p,tp,pp))))
