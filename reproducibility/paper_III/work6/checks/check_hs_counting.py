"""Finite-rank multiplicity/check of Dirac vs Weyl prefactors using a toy F.
The toy F is explicitly NOT the Schwarzschild transition function.
No inference about exact Schwarzschild uniform angular tails is made.
"""
from pathlib import Path
import csv,math,json
import numpy as np
from scipy.integrate import quad
BASE=Path(__file__).resolve().parent

def ftoy(x):
    return .5/(1+x**8)

def sum_at_k(k,window=None):
    upper=int(math.ceil((20 if window is None else window[1])*k**.25))+2
    lambdas=np.arange(1,upper+1,dtype=float)
    xs=lambdas/k**.25
    mask=np.ones_like(xs,dtype=bool) if window is None else (xs>=window[0])&(xs<=window[1])
    return float(np.sum(4*lambdas[mask]*ftoy(xs[mask])))

I, _=quad(lambda x:x*ftoy(x),0,np.inf,epsabs=1e-13)
I_exact=math.pi/(8*math.sqrt(2))
assert abs(I-I_exact)<1e-12
Ia,_=quad(lambda x:x*ftoy(x),.25,.75,epsabs=1e-13)
L=2*math.pi
rows=[]
for K in (100,1000,10000):
    sk_win=2*sum(sum_at_k(k,(.25,.75)) for k in range(1,K+1))
    sk_all=2*sum(sum_at_k(k) for k in range(1,K+1))
    row=dict(K=K,normalized_window=sk_win/(L*K**1.5),expected_window=8*Ia/(3*math.pi),normalized_full=sk_all/(L*K**1.5),expected_full=8*I/(3*math.pi))
    rows.append(row);print(row)
# assert sufficient convergence in the large-K toy benchmark
assert abs(rows[-1]['normalized_full']/rows[-1]['expected_full']-1)<.03
assert abs(rows[-1]['normalized_window']/rows[-1]['expected_window']-1)<.04
with (BASE/'toy_HS_counting.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
(BASE/'counting_validation.json').write_text(json.dumps({'I_full_toy':I,'I_full_toy_closed':I_exact,'I_window_toy':Ia,'rows':rows,'note':'Toy benchmark validates count of both k signs, Dirac two polarizations and 2lambda angular multiplicity. Does not certify physical transition F.'},indent=2))
print('Checks OK')
