#!/usr/bin/env python3
"""Fast independent verification subset for Paper II v1.2.0. Original audit scripts unchanged."""
import importlib.util
from pathlib import Path
import numpy as np

source = Path(__file__).resolve().parent / "audit_num.py"
spec = importlib.util.spec_from_file_location("paper_II_independent_audit_num", source)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

for l, K, target in [(1,1e4,1.00864), (2,1e5,1.04792), (3,1e7,2.16185)]:
    r0 = mod.r_of_eta(200/K)
    a,b,wr,_ = mod.integrate_jost(K,l*(l+1),r0)
    al,be,_ = mod.endpoint_coeffs(a,b,K,l*(l+1))
    val = abs(be)**2
    assert abs(val-target)<7e-6, (l,K,val,target)
    assert abs(abs(al)**2-abs(be)**2-1)<3e-10, (l,K)
    print(f"PASS l={l} K={K:.0f} n={val:.9f}")

for lh,kap,target in [(0.5,10.0,-0.44016),(0.8,40.0,+0.36638)]:
    L2 = lh*kap**2
    K = np.sqrt(kap**2-L2)
    a,b,wr,_ = mod.integrate_jost(K,L2,1-1e-9)
    al,be,_ = mod.endpoint_coeffs(a,b,K,L2)
    bench = 1/(36*np.pi**2*lh**2)
    rel = 100*(abs(be)**2*kap**4/bench-1)
    assert abs(rel-target)<0.003, (lh,kap,rel,target)
    print(f"PASS lhat2={lh} kappa={kap:.0f} directional_deviation_pct={rel:+.4f}")
print("PAPER II INDEPENDENT FAST CHECK PASS (not full numerical audit)")
