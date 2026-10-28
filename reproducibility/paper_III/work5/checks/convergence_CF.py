"""Convergence campaign for universal-model C_F (not Schwarzschild full HS coefficient).
C_F=integral_0^infty u F(u^(4/3))du. Numerical cutoff at Lambda_c;
leading analytic tail is (9/4096)/6 * Lambda_c^(-9/2).
A tail estimate, not a rigorous enclosure, since O(Lambda^-7) coefficient unknown.
"""
import numpy as np, csv, time, json
from scipy.special import roots_legendre
from universal_precision import transition
from pathlib import Path
DIR=Path(__file__).parent
rows=[]
cutoff=8.
U=cutoff**.75
TAIL=(9/4096)/6 * cutoff**(-4.5)
with (DIR/'CF_quadrature_convergence.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['cutoff_Lambda','U','N_gauss','C_num','tail_leading','C_estimate','max_abs_norm_error','sum_nfev','wall_seconds']);w.writeheader()
    for N in [12,20,32,48,64,96]:
        start=time.monotonic(); xi,wi=roots_legendre(N)
        u=.5*U*(xi+1);w_u=.5*U*wi
        sf=0.; maxnorm=0.; ne=0
        for i,(uv,ww) in enumerate(zip(u,w_u)):
            Lm=float(uv**(4/3))
            tr=transition(Lm,rtol=3e-12,order=2)
            sf += float(uv*ww*tr['F']); maxnorm=max(maxnorm,abs(tr['norm_error']));ne+=tr['nfev']
        row={'cutoff_Lambda':cutoff,'U':U,'N_gauss':N,'C_num':sf,'tail_leading':TAIL,'C_estimate':sf+TAIL,'max_abs_norm_error':maxnorm,'sum_nfev':ne,'wall_seconds':time.monotonic()-start}
        print(row,flush=True);w.writerow(row);f.flush();rows.append(row)
(DIR/'CF_quadrature_summary.json').write_text(json.dumps(rows,indent=2))
