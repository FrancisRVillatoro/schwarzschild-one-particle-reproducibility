import numpy as np, math
from num_checks import *
D = np.pi**1.5*math.gamma(0.25)/(24*math.gamma(0.75)**3)
print('s wave, phase check: beta*sqrt(K) vs D(1+i)=%.4f(1+i)'%D)
for K in (32, 128, 512):
    al, be = endpoint_coeffs(K, 0); print('  K=%4d  beta*sqrt(K)='%K, np.round(be*np.sqrt(K),4), ' K|beta|^2=%.6f'%(K*abs(be)**2))
print('validation of LG start against full integration, l=1, K=1000:')
a1,b1 = endpoint_coeffs(1000,1); a2,b2 = endpoint_coeffs_LG(1000,1)
print('  full alpha,beta=',np.round(a1,6),np.round(b1,6),'  LG=',np.round(a2,6),np.round(b2,6))
print('\n l        K     |b|^2 num   |b|^2 asym   rel.dev     alpha (num)           alpha (thm)          beta (num)            beta (thm)')
rows=[]
for l in (1,2,3):
    for K in (1e2,1e3,1e4,1e5,1e6,1e7):
        al,be = endpoint_coeffs(K,l) if K<1000 else endpoint_coeffs_LG(K,l)
        Lam = K/(4*np.exp(3*EG)*(l*(l+1))**2); lt = np.log(Lam)/(2*np.pi)
        ap, bp = 1.25-1j*lt, -0.75-1j*lt
        nn, na = abs(be)**2, abs(bp)**2
        rows.append((l,K,nn,na))
        print(f'{l:2d} {K:9.0e}  {nn:10.5f}  {na:10.5f}  {(nn-na)/na:+.2e}   {al.real:+.4f}{al.imag:+.4f}i   {ap.real:+.4f}{ap.imag:+.4f}i   {be.real:+.4f}{be.imag:+.4f}i   {bp.real:+.4f}{bp.imag:+.4f}i   norm-1={abs(al)**2-abs(be)**2-1:+.0e}')
np.save('part1_rows.npy', np.array(rows))
