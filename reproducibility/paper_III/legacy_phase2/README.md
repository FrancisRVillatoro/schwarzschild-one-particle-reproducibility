# Historical phase-0/phase-2 fermionic diagnostics

These supplemental Paper III diagnostics preserve the historical source files.
The Work1–Work8 sources provide the production checks.

`dirac_phase0.py` supplies H, projections and evolution used by
`dirac_phase2_numerics.py` and `dirac_phase2_packet.py`.
`dirac_stress_identity_check.py` is independent of those imports.

## Conventions and limits

The legacy evolution uses i dchi/dsigma=H chi, opposite to the physical future
orientation in the current paper. For the real reduced matrices, simultaneous
conjugation of coefficients and evolution preserves the mode occupations. It
is not a replacement for the current physical two-point kernel. Coherent packet
phases cannot be transferred unchanged to the physical convention.

The two-frequency packet is a finite-mode auxiliary diagnostic, not a normalized
physical-continuum Unruh excitation. Use work3/work4 for continuum normalization,
CAR, Ward identities and all stress components; do not infer positivity of the
stress or a bound on the whole tensor from this legacy packet.

The supplied counts and packet output files are retained. The phase-0 standalone
output and stress-identity output are new reruns recorded during integration.
The counts reproduce all printed occupations; the last residual at one point
can change by roundoff (4.48e-15 versus 4.45e-15). The packet output reproduces
byte for byte in the recorded environment.

The refined [Work5 scripts and reference outputs](../work5/checks/) support the
small-Lambda results in the manuscript.
