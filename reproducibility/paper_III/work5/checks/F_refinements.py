"""Refinements of Jost initialization, integration tolerance and finite starting radius.
Floating-point convergence diagnostics, not interval certification.
"""
from universal_precision import transition
import csv, json
from pathlib import Path
D=Path(__file__).parent
vals=[0.0001,0.001,0.01,0.02,0.05,0.1,0.2,0.3,0.5,0.7,1,1.3,1.6,2,3,4,4.4,4.6,5,6,8,10,12]
rows=[]
with (D/'F_refinements.csv').open('w',newline='') as f:
    names=['Lambda','rho0_base','F_order0','F_order2','F_refined','difference_start','difference_order','nfev_refined','norm_defect_refined']
    w=csv.DictWriter(f,fieldnames=names);w.writeheader()
    for L in vals:
        # nominal rho0 yields ~1024 units in lambda*rho0^2, or 20 at high lambda
        rho=max(20.,32/L**.5)
        raw=transition(L,rho0=rho,rtol=3e-12,order=0)
        base=transition(L,rho0=rho,rtol=3e-12,order=2)
        fine=transition(L,rho0=1.4*rho,rtol=1e-12,order=2)
        row=dict(Lambda=L,rho0_base=rho,F_order0=raw['F'],F_order2=base['F'],F_refined=fine['F'],difference_start=abs(fine['F']-base['F']),difference_order=abs(base['F']-raw['F']),nfev_refined=fine['nfev'],norm_defect_refined=fine['norm_error'])
        w.writerow(row); f.flush();rows.append(row);print(row,flush=True)
(D/'F_refinements.json').write_text(json.dumps(rows,indent=2))
