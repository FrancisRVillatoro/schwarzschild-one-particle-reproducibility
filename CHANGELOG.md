# Changelog

## v1.3.0 — 2026-11-10

- Added the Paper II supplementary scalar scripts for three s-wave refinements, 48/96-node packet refinements, the normalized Table 2 and 20 independent Ward-identity cases.
- Exported the selected Table 1 and Table 2 rows with their numerical settings and source hashes.
- Limited the independent scalar stdout tolerance to absolute differences of at most `1e-13` in the B1/B3 norm residuals. The other fields remain exact after elapsed times are removed.
- Bounded the infrared Jost norm defect by `1e-6`, with exact frequencies and the standard tolerances on the scaled amplitude.
- Checked the derived `additive_error` and `bound_ratio` in the Paper III terminal-bound rows by reconstruction with a roundoff budget.
- Added historical fermionic phase-0 and phase-2 diagnostics with their auxiliary convention documented.
- Verified root and nested manifest coverage and freshly generated fermionic tables against the frozen references. Full validation passed in WSL.
- Preserved the v1.1.0 scientific sources and their published reference outputs.

## v1.2.0 — 2026-10-28

- Aligned Paper II metadata and updated its manuscript cross-references.
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
