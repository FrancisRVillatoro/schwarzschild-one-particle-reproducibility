"""Numerical diagnostic for the window HS coefficient (universal F, no Schwarzschild UV-tail claim)."""
from pathlib import Path
import numpy as np, csv
from scipy.special import roots_legendre
from universal_precision import transition
out=Path(__file__).resolve().parent/'window_integrals.csv'
a,b=0.25,0.75
with out.open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['a','b','N','I_window','8_over_3pi_times_I','max_norm_defect','nfev'])
    for N in (8,16,24,32):
        z,weights=roots_legendre(N)
        u=a+(b-a)*(z+1)/2
        ww=weights*(b-a)/2
        vals=[]; norm=[]; nfev=0
        for v in u:
            r=transition(float(v**(4/3)),rtol=3e-12,order=2)
            vals.append(r['F']);norm.append(abs(r['norm_error']));nfev+=r['nfev']
        i=float(np.dot(ww,u*np.array(vals)))
        row=[a,b,N,f'{i:.15g}',f'{i*8/(3*np.pi):.15g}',f'{max(norm):.4g}',nfev]
        print(row,flush=True);w.writerow(row);f.flush()
