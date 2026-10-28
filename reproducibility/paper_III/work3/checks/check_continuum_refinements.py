#!/usr/bin/env python3
"""Frequency quadrature and ODE refinement of the continuum diagnostic (not Unruh)."""
import json,csv
from pathlib import Path
from check_continuum_stress import continuum_check
rows=[]
for N,rtol,atol in [(48,2e-12,2e-14),(96,2e-12,2e-14),(144,2e-12,2e-14),(96,2e-13,2e-15)]:
 o=continuum_check(N,rtol,atol,save=False)
 rows.append({'N_kappa':N,'rtol':rtol,'atol':atol,'rho_limit':o['rho_leading_coefficient'],'shear_limit':o['T12_leading_coefficient'],'scaled_rho_r_1e_6':o['rows'][-1]['scaled_rho'],'scaled_shear_r_1e_6':o['rows'][-1]['scaled_shear'],'max_norm_error':o['max_norm_error']})
p=Path(__file__).with_name('continuum_refinements.csv')
with p.open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps(rows,indent=2))
