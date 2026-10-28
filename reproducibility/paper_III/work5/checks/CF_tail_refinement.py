"""Resolve C_F beyond Lambda=8 by direct quadrature of the universal ODE.
Asymptotic tail is estimated, not an interval-rigorous bound.
"""
from pathlib import Path
import numpy as np, csv,json,time
from scipy.special import roots_legendre
from universal_precision import transition
DIR=Path(__file__).parent

with (DIR/'CF_quadrature_convergence.csv').open() as f:
    ref=list(csv.DictReader(f))
q=float(next(row['C_num'] for row in ref if int(row['N_gauss']) == 96)) # computed 0<Lambda<8
rows=[]
for lo, hi,N in [(8.,10.,16),(10.,12.,16),(12.,16.,20)]:
    t=time.monotonic(); ul=lo**.75;uh=hi**.75
    xi,wi=roots_legendre(N);uv=ul+(uh-ul)*(xi+1)/2; wv=(uh-ul)*wi/2
    piece=0.; nfev=0
    for u,w in zip(uv,wv):
        tr=transition(float(u**(4/3)),rtol=3e-12,order=2)
        piece+=float(w*u*tr['F']); nfev+=tr['nfev']
    q+=piece
    cutoff=hi
    tail=(9/4096)/6 * cutoff**(-4.5)
    row=dict(Lambda_interval_lo=lo,Lambda_interval_hi=hi,N_gauss=N,interval_integral=piece,C_numerical_to_cutoff=q,leading_tail=tail,C_estimate=q+tail,wall_time=time.monotonic()-t,nfev=nfev)
    rows.append(row);print(row,flush=True)
    with (DIR/'CF_cutoff_convergence.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(row));writer.writeheader();writer.writerows(rows)
(DIR/'CF_cutoff_summary.json').write_text(json.dumps(rows,indent=2))
