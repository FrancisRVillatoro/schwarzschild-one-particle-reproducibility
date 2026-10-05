import numpy as np, sys
sys.argv=['x','none']
exec(open('convergence_study.py').read().split("which=sys.argv[1]")[0])
from numpy.polynomial.legendre import leggauss
def B(L,Nu,target,rtol=1e-11):
    xu,wu=leggauss(Nu); us=(xu+1)/2; wu=wu/2
    return sum(wU*(3*u*u-1)*exact_A(L,float(u),target,rtol) for u,wU in zip(us,wu))+1/(63*L**3)
d7=4/2835
Ls=np.arange(8.,16.01,1.)
y=np.array([B(L,48,360.)*L**7-d7 for L in Ls])
np.save('tail_y.npy',np.vstack([Ls,y]))
def tail_from(d,cut=16.):   # d = {k: d_k}, radial integrand L^4 B /(4 pi^2)
    return sum(dk/((k-5)*cut**(k-5)) for k,dk in d.items())/(4*np.pi**2)
res={}
for name,mask,ks in [('2t[8,16]',Ls>=8,(9,11)),('3t[8,16]',Ls>=8,(9,11,13)),('2t[10,16]',Ls>=10,(9,11)),('3t[10,16]',Ls>=10,(9,11,13))]:
    M=np.column_stack([Ls[mask]**(7-k) for k in ks]); c=np.linalg.lstsq(M,y[mask],rcond=None)[0]
    d={7:d7}; d.update({k:ck for k,ck in zip(ks,c)})
    res[name]=tail_from(d); print(name,'coeffs',dict((k,float(v)) for k,v in d.items()),' tail16=',res[name])
print('baseline tail16 (archive method) = 7.055499424007197e-08 ; d7-only tail =',tail_from({7:d7}))
