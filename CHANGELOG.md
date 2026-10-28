# Changelog

## v1.2.0 — 2026-10-28

- Aligned Paper II metadata with its Annals of Physics submitted TEX and updated its MANIFEST cross-references.
- Added `reproducibility/paper_II/independent_audit/`: independent symbolic and numerical audit scripts, author-supplied outputs, SHA checks, and quick regression test.
- Preserved the Paper-I and Paper-II production scripts from v1.1.0 without changes.

- Added `reproducibility/paper_III` for the Weyl/Dirac study.
- Included exact algebraic CAR, spin-harmonic and symbol checks, plus massive-mode and stress tests.
- Added reproducible `F(Lambda)` refinement CSVs that show nonmonotonicity, and `C_F` cutoff/quadrature studies.
- Kept full model `C_F` distinct from exact Schwarzschild fixed-window Hilbert--Schmidt bounds.
- Added fast and full Paper-III runners and independent provenance/hash manifests.
- Preserved all v1.1.0 Paper-I/Paper-II production scripts and reference results.

## v1.1.0 — 2026-10-05

- Reorganized reproducibility material into separate `paper_I` and `paper_II` trees.
- Added Paper-I convergence studies for `C_dyn` and `Delta_dyn` and revised uncertainty estimates.
- Replaced the earlier one-panel estimate of the subleading RSET ultraviolet tail by a direct fit including `lambda^-9`, `lambda^-11` and stability variants.
- Added exact symbolic adiabatic checks used in the revised Paper I.
- Added equal-time Hadamard point-splitting calculation on the Kasner geometry.
- Added Paper-II symbolic and numerical checks for endpoint-Hankel asymptotics, fixed-angular towers, generic rays, exceptional-slice normalization and one-particle stress.
- Retained selected v1.0.0 outputs explicitly marked as legacy for auditability.

## v1.0.0 — 2026-08-19

Initial public reproducibility release.
