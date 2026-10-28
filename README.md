# Schwarzschild one-particle and Kasner fixed-point reproducibility

Numerical and symbolic reproducibility package supporting three related manuscripts on quantum fields near the Schwarzschild spacelike singularity:

1. Francisco R. Villatoro, *A universal Kasner scaling limit for Hadamard states in the Schwarzschild interior* (2026).
2. Francisco R. Villatoro, *Modewise freezing, terminal particle structures and local stress in the Schwarzschild interior* (2026).
3. Francisco R. Villatoro, *Weyl and Dirac fields near the Schwarzschild singularity: terminal structure, universal spin rotation and local stress* (2026).

The repository contains reproducibility material only: scientific Python code, symbolic checks, compact reference outputs, and documentation. It does not contain manuscript source or PDFs.

## Version 1.1.0

Version 1.1.0 reorganizes the package into independent `paper_I/` and `paper_II/` reproducibility trees and adds the calculations used in the revised manuscripts.

Major additions relative to v1.0.0 include:

- convergence studies and uncertainty estimates for `C_dyn` and `Delta_dyn`;
- direct ultraviolet fitting of the subleading `lambda^-9` and higher terms for the Kasner RSET;
- exact symbolic adiabatic checks used in Paper I;
- equal-time Hadamard point splitting on the Kasner geometry and the finite constants `1/(240*pi^2)` and `-523/(498960*pi^2)`;
- numerical and symbolic checks for the Paper II endpoint-Hankel asymptotics, fixed-angular towers, generic rays, the exceptional maximal-volume slice, and one-particle stress asymptotics.

The archived v1.0.0 release remains immutable and can be used to reproduce the earlier manuscript state.


## Paper II — independent Annals of Physics audit (v1.2.0)

The v1.1.0 scalar production scripts remain unchanged. Version v1.2.0 adds
`reproducibility/paper_II/independent_audit/` with separate symbolic and
numerical scripts, submitted reference outputs, SHA-256 checks and a short
independent fast reproduction. Run
`python reproducibility/paper_II/independent_audit/check_quick_independent.py`
for fast verification; the full `audit_num.py` includes an expensive
horizon-start test at K=10^4 and special-function quadratures and can take
substantially longer. The versioned archive does not modify the Paper II
manuscript or v1.1.0. Consult `reproducibility/MANIFEST.md` for scope.

## Paper III — Weyl and Dirac fields (v1.2.0)

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
- `reproducibility/paper_II/` — original scalar production checks plus a v1.2.0 independent symbolic/numerical audit.
- `reproducibility/paper_III/` — fermionic CAR/Hadamard, spin-harmonic, universal rotation, Hilbert--Schmidt and continuum-stress checks.
- `reproducibility/MANIFEST.md` — detailed mapping from manuscript results to scripts and frozen outputs.
- `run_quick_checks.sh` — fast symbolic checks for both papers.
- `run_full_reproduction.sh` — principal numerical and symbolic calculations for both papers.
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

The Jiang--Jiang Schwarzschild coefficients are used only as external benchmarks. No manuscript source or PDF is stored in this repository. The manuscripts are scientifically distinct but share this reproducibility archive.

## License

The code and associated reproducibility files are released under the MIT License.

## Funding

The earlier scalar reproducibility project acknowledges PPRO-B4-2026-004, Universidad de Málaga. Funding applicable to the fermionic study must be confirmed by the author before use in the journal manuscript.
