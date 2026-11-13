# Paper III — fermionic reproducibility

Supporting *Weyl and Dirac fields near the Schwarzschild singularity: terminal structure, universal spin rotation and local stress* (2026), separate from the two scalar papers.

This folder contains the exact scripts and frozen reference outputs collected from successive Work1--Work8 audits. **These outputs are inherited from the earlier work bundles.** To claim a fresh reproduction, run the scripts using the runners below and retain the new logs/reports. The checks test algebra, numerical approximations, and explicit examples, not the infinite-dimensional theorems.

## Reproduction

The environment used for the preceding runs contained Python 3.11--3.13, NumPy 2.3.5, SciPy 1.17.0, SymPy 1.14.0 and mpmath 1.3.0. Install the root `requirements.txt` in a fresh virtual environment (Python 3.11+).

From the repository root (runs in an isolated copy, leaving the tracked CSV/JSON unchanged):

```bash
bash run_quick_checks_paper_III.sh
bash run_full_reproduction_paper_III.sh
```

The fast suite runs the time-orientation, symbol-bookkeeping, spin/algebra, CAR and massive operator checks, plus a freshly integrated nonmonotonic `F` window. The 14-stage full suite also regenerates the invariant angular kernels, two continuum-packet campaigns, Ward/CAR diagnostics, the detailed `F` and `C_F` convergence studies, the window coefficient and the Hilbert--Schmidt toy multiplicities. **Full mode is CPU-intensive.**

For each row of `work4/checks/terminal_bounds.csv`, the full runner checks that
the additive density error lies within the explicit bound and reconstructs
`additive_error` and `bound_ratio` from the reported quantities. It allows
roundoff before subtracting the two density terms and propagates that budget
through the ratio. The report records the row-level budgets and comparisons.

## Scientific mapping

- `work1/`: tetrad time orientation, Clifford/CAR signs and mode evolution.
- `work2/checks/`: algebraic superadiabatic identities; the all-order argument is in the paper.
- `work3/checks/`: spin spherical kernel and continuum stress diagnostics.
- `work4/checks/`: Ward identities, restricted CAR and continuum conservation.
- `work5/checks/`: `F(Lambda)`, superadiabatic starting data, `C_F` quadratures.
- `work6/checks/`: Hilbert--Schmidt **fixed-window** integral and toy multiplicities.
- `work7/`: massive operator-block diagnostics.
- `work8/`: exact CAR purity/restriction counterexample and rank-one excitations.

### Critical distinctions

1. `C_F ~= 0.07888494` is a constant **of the universal spin-rotation model**. It is **not** claimed as the complete angular Schwarzschild Hilbert--Schmidt coefficient. The Schwarzschild theorem establishes fixed-window equivalence and a global lower bound.
2. `F(Lambda)` is **not globally monotone**; the frozen file `work5/checks/F_refinements.csv` includes Lambda=4.4,4.6,5.0.
3. The continuum stress has leading comoving density of order `R^-5/2` only when its packet coefficient is nonzero; mixed spatial spin-connection components may be `R^-3`.
4. The continuum stress packet tests use explicitly prescribed Cauchy/terminal data. They are **not a numerical computation of the global Unruh transmission matrix**.
5. The absolute renormalized Dirac vacuum stress and its transfer from a distributional scaling limit are not calculated.
6. The global-state-to-interior CAR restriction is generally mixed and contractive; the transmitted interior packet is **not renormalized** to norm one.

## Release and version

This guide describes [v1.3.0](../../CITATION.cff). Cite the release tag and its
version-specific Zenodo DOI once the archive is available. The local
`SHA256SUMS` covers the Paper III tree.
