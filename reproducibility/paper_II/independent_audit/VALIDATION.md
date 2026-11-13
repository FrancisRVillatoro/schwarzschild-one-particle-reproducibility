# Historical independent-audit execution

In the earlier independent-audit execution, both supplied scripts were run
unmodified, including all parts of audit_num.py
(no SKIP_A). audit_sym stdout is byte-identical to the supplied output.
The complete audit_num stdout is identical after removing only elapsed-second
suffixes; no scientific number or residual was masked. Original scripts and output
files are kept byte for byte. Environment: Python 3.13.5, NumPy 2.3.5,
SciPy 1.17.0, SymPy 1.14.0, mpmath 1.3.0. The
[audit guide](README.md) describes the v1.3.0 comparison policy.

Matching output is not certification of every printed benchmark. In particular,
the supplied epsilon=0.02 Abel quadrature and its Richardson extrapolation are not
converged to 2/pi^2. The printed extrapolation is 0.20210119-0.00060915263i,
whereas 2/pi^2=0.20264237. This historical diagnostic is not promoted to a high-
accuracy proof or a validated constant. No code/output was silently repaired.
The exact identities and production symbolic audit remain separately identified.
