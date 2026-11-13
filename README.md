# Schwarzschild one-particle and Kasner fixed-point reproducibility

Numerical and symbolic reproducibility package supporting three related manuscripts on quantum fields near the Schwarzschild spacelike singularity:

1. Francisco R. Villatoro, *A universal Kasner scaling limit for Hadamard states in the Schwarzschild interior* (2026).
2. Francisco R. Villatoro, *Particle creation, terminal states and local stress of a scalar field in the Schwarzschild interior* (2026).
3. Francisco R. Villatoro, *Weyl and Dirac fields near the Schwarzschild singularity: terminal structure, universal spin rotation and local stress* (2026).

The repository contains reproducibility material only: scientific Python code, symbolic checks, compact reference outputs, and documentation. It does not contain manuscript source or PDFs.

## Version 1.3.0 — 2026-11-10

Version 1.3.0 adds the Paper II supplementary scalar calculations, refined
publication tables and Ward-identity checks to the independent scalar audit and
fermionic suite introduced in v1.2.0. It includes historical phase-0 and phase-2
fermionic diagnostics alongside the Work1–Work8 checks. The v1.1.0 scientific
sources and reference outputs are preserved byte for byte.

The joint validator allows an absolute difference of `1e-13` in the printed
B1/B3 scalar norm residuals. It checks the infrared Jost norm defect against an
absolute bound of `1e-6` and compares the scaled amplitude at the standard
CSV tolerances. The fermionic runner reconstructs the two derived error columns
in `terminal_bounds.csv` and accounts for roundoff in their comparison. See the
[reproducibility manifest](reproducibility/MANIFEST.md) for the result mapping
and exact refinement settings.

Full validation passed in WSL, including the independent scalar audits, the
supplementary scalar table and Ward checks, and all 14 fermionic stages. The
report records the checks and comparisons and identifies any saved audit output
reused during a resumed run. Run the joint checks from the repository root:

```bash
python3 tools/validate_v120.py --mode quick --report /tmp/schwarzschild_v130_quick.json
python3 tools/validate_v120.py --mode full --report /tmp/schwarzschild_v130_full.json
```

Use Python 3.11 or newer with the root requirements installed. The validator
checks the release manifests and runs the scripts in temporary copies. Full
mode includes the complete numerical audit with `SKIP_A` unset.

Cite the GitHub release tag [v1.3.0](https://github.com/FrancisRVillatoro/schwarzschild-one-particle-reproducibility/releases/tag/v1.3.0)
and the version-specific Zenodo DOI once the archive is available.

## Version 1.2.0 — 2026-10-28

The [v1.2.0 release](https://github.com/FrancisRVillatoro/schwarzschild-one-particle-reproducibility/releases/tag/v1.2.0)
introduced the Paper II independent symbolic and numerical audit and the Paper
III Weyl/Dirac reproducibility tree, including the fast and full fermionic
runners.

## Version 1.1.0

Version 1.1.0 reorganizes the package into independent `paper_I/` and `paper_II/` reproducibility trees and adds the calculations used in the revised manuscripts.

Major additions relative to v1.0.0 include:

- convergence studies and uncertainty estimates for `C_dyn` and `Delta_dyn`;
- direct ultraviolet fitting of the subleading `lambda^-9` and higher terms for the Kasner RSET;
- exact symbolic adiabatic checks used in Paper I;
- equal-time Hadamard point splitting on the Kasner geometry and the finite constants `1/(240*pi^2)` and `-523/(498960*pi^2)`;
- numerical and symbolic checks for the Paper II endpoint-Hankel asymptotics, fixed-angular towers, generic rays, the exceptional maximal-volume slice, and one-particle stress asymptotics.

The archived v1.0.0 release remains immutable and can be used to reproduce the earlier manuscript state.


## Paper III — Weyl and Dirac fields

Fermionic reproducibility source and frozen numerical checks have been added in
`reproducibility/paper_III/`. They support the separate manuscript
*Weyl and Dirac fields near the Schwarzschild singularity: terminal structure, universal spin rotation and local stress*.
The archive includes unitary terminal evolution and Clifford/CAR checks; invariant
spin-harmonic kernels; angular-kernel convergence; the universal two-level
transition function `F(Lambda)` including nonmonotonicity and sudden/adiabatic
limits; finite-window Hilbert--Schmidt counting; the massive transition block;
and physical-continuum stress, box-limit and Ward-identity checks.

The numerically evaluated constant `C_F ~= 0.07888494` belongs to the
**universal two-level model**, not to the complete Schwarzschild angular sum.
The fermionic scripts do **not** compute the vacuum Dirac RSET or an actual
exterior Unruh transmission matrix. See `reproducibility/paper_III/README.md`.

For a fast standalone check of Paper III, use
`bash run_quick_checks_paper_III.sh`. For the full numerical campaign, use
`bash run_full_reproduction_paper_III.sh`. The original Paper-I/II commands
remain unchanged from v1.1.0.

## Repository layout

- `reproducibility/paper_I/` — Kasner scaling fixed-point numerical/symbolic calculations and convergence audits.
- `reproducibility/paper_II/` — particle-number, endpoint-Hankel, ultraviolet-obstruction and one-particle-stress checks.
- `reproducibility/paper_III/` — fermionic CAR/Hadamard, spin-harmonic, universal rotation, Hilbert--Schmidt and continuum-stress checks.
- `reproducibility/MANIFEST.md` — detailed mapping from manuscript results to scripts and frozen outputs.
- `run_quick_checks.sh` — fast symbolic checks for Papers I and II.
- `run_full_reproduction.sh` — principal numerical and symbolic calculations for Papers I and II.
- `tools/validate_v120.py` — joint validation for the Paper II supplementary calculations and the Paper III suite.
- `requirements.txt` — Python dependencies.
- `CITATION.cff` — citation metadata used by GitHub/Zenodo.
- `SHA256SUMS` — release integrity manifest.

## Environment

Install the pinned dependencies with

```bash
python -m pip install -r requirements.txt
```

The scripts require NumPy, SciPy, pandas, SymPy and mpmath.

## Quick checks

From the repository root:

```bash
bash run_quick_checks.sh
```

These rederive the main exact symbolic identities and the Kasner point-splitting constants for Paper I and the principal symbolic identities of Paper II.

## Full numerical reproduction

From the repository root:

```bash
bash run_full_reproduction.sh
```

The convergence campaign is documented by the frozen files `reproducibility/paper_I/conv_C.txt`, `conv_D.txt`, `convergence_results.txt` and `tail_y.npy`. See `reproducibility/MANIFEST.md` for the exact commands, quoted uncertainties and run times.

The final Paper-I numerical values archived in v1.1.0 are

```text
C_dyn     = -2.626175(4)e-4
C_K       =  1.595541(4)e-4
Delta_dyn =  7.06092(2)e-5
Delta_K   = -3.55936(2)e-5
```

## Production versus legacy files

Files explicitly marked as legacy in `reproducibility/MANIFEST.md` are retained only to document the v1.0.0 calculation. The revised production result for Paper I is generated by `reproducibility/paper_I/kasner_complete_rset.py`; it uses the convergence study and the direct ultraviolet tail fit that includes the `lambda^-11` contribution.

## Data availability and scope

The Jiang--Jiang Schwarzschild coefficients are used only as external benchmarks. No manuscript source or PDF is stored in this repository. The three manuscripts are scientifically distinct but share this reproducibility archive.

## License

The code and associated reproducibility files are released under the MIT License.

## Funding

The scalar reproducibility project acknowledges PPRO-B4-2026-004, Universidad de Málaga.
