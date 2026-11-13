# Reproducibility package for the three Schwarzschild-interior papers

Version 1.3.0. Paper I entries retain their v1.1.0 mapping. Use Python 3.11 or
newer with the pinned dependencies in [requirements.txt](../requirements.txt).
Paper I uses G = c = hbar = 1 and Kasner coordinates with y = 1. Paper II uses
2M = 1 unless stated otherwise. The adjacent output files retain the v1.1.0
scalar records. The supplementary scripts also write the CSV and JSON files
listed below. The joint validator runs in disposable copies.

## paper_I  (A universal Kasner scaling limit for Hadamard states in the Schwarzschild interior)

| Manuscript item | Script | Output file | Run time |
|---|---|---|---|
| Eq. (5.5) C_dyn baseline mode integral; UV tail of the bracket (c7 exact, c9, c11 fitted) | `reproduce_kasner_phi2.py` (uses nothing else) | `kasner_phi2_reproduction.csv` (archive v1.0.0 run) | 45 s |
| Table 2 (convergence of C_dyn and Delta_dyn), final values C_dyn = -2.626175(4)e-4, Delta_dyn = 7.06092(2)e-5 | `convergence_study.py C <variant...>` and `D <variant...>`; `tailfit.py` (direct fit of d9, d11 on [8,16]) | `conv_C.txt`, `conv_D.txt`, `tail_y.npy`, `convergence_results.txt` | 40-150 s per variant |
| Eqs. (5.6)-(5.7), (5.9), (5.11)-(5.12), Table 1 (all Kasner and Schwarzschild coefficients, conformal values, null-projection thresholds, Jiang-Jiang comparison) | `kasner_complete_rset.py` (needs `kasner_rset_core.py`, `convergence_study.py`, `reproduce_kasner_phi2.py` in the same directory) | `kasner_complete_rset_output.txt` | 150 s |
| Appendix A: identity (5.6), residue 1/(240 pi^2), J_mu of Eq. (A.5), residue -523/(498960 pi^2), vanishing massless fourth-order average, sixth-order tails -1/2835 and 4/2835, polynomials B_2 and C_4; Section 5: coupling shifts (2,-10/3,-16/3) and (9,-15,-24), Kretschmann 64/27 and 48 M^2/r^6 | `kasner_adiabatic_symbolics.py` | `kasner_adiabatic_symbolics_output.txt` | 3 min |
| Appendix C and Proposition on Hadamard renormalization on Kasner: sigma, U, V_0, V_1 from their defining equations, equal-time polynomials (C.3), finite parts, constants 1/(240 pi^2) and 523/(498960 pi^2), matching of singular terms | `kasner_point_splitting.py` | `kasner_point_splitting_output.txt` | 2 min |
| Section 5.2 and Appendix A, mass regulator: exact massive modes versus massless, infrared layer | `check_kasner_massive_IR.py` (archive v1.0.0) | `kasner_massive_IR_layer.csv`, `kasner_massive_IR_renormalized_layer.csv` | 10 min |
| Legacy checks from archive v1.0.0 (kept for the record) | `check_kasner_rset_symbolics.py`, `check_kasner_uv_tails.py` (the latter only prints c7, d7; their derivation is in `kasner_adiabatic_symbolics.py`); `kasner_delta_dynamic_pieces.csv`, `kasner_complete_rset_summary_v1_0_0.csv` (values superseded by `kasner_complete_rset_output.txt`: the earlier tail of Delta_dyn omitted the lambda^-11 term) | | |

## paper_II — Particle creation, terminal states and local stress of a scalar field in the Schwarzschild interior

Equation and table labels follow the Paper II manuscript. All paths in this
block are relative to `reproducibility/paper_II/`. The five v1.1.0 production
scripts and their published outputs remain the scalar baseline.

| Manuscript item | Source | Reference output / scope |
|---|---|---|
| Table 1 upper block, massless s wave; original default run | `num_checks.py` (unchanged) | `num_checks_output.txt`; use the refined supplemental rows for revised printed digits |
| Proposition 9, Eqs. (6.3)–(6.4), phase; Table 1 lower block and Theorem 10, Eqs. (6.8)–(6.10) | `run_part1.py` (unchanged; imports `num_checks.py`) | `run_part1_output.txt` |
| Table 2 original raw horizon diagnostic | `run_part2.py` (unchanged; imports `num_checks.py`) | `run_part2_output.txt`; no normalization silently added to this historical output |
| Eq. (6.6), directional first-Born benchmark | `generic_ray_check.py` (unchanged) | `generic_ray_check_output.txt`; benchmark, not a proof of a uniform ray remainder |
| Eq. (2.3), pump; Proposition 2, velocity map and r=3/4; Theorem 10 algebra; Remark 18, coupling | `paper2_symbolic_checks.py` (unchanged) | `paper2_symbolic_checks_output.txt` |
| Independent symbolic re-derivation: pump (2.3), rescaling (3.5), Lemma 5, series (6.2), geodesics (8.9), Frobenius Appendix F, counts (3.2),(3.4),(3.6), Remark 18 | `independent_audit/audit_sym.py` | `independent_audit/audit_sym_output.txt`; reproduced exactly |
| Independent integrator: Table 1, Proposition 9 integrals/constants, Theorem 10 towers, Eq. (6.6) | `independent_audit/audit_num.py` | `independent_audit/audit_num_output.txt`; exact comparison apart from elapsed times and the B1/B3 norm-residual tolerance |
| Table 1 upper block: three endpoint/tolerance levels; Table 2: Z-normalized profile and 48/96 nodes; 29 additional checks | `aop_supplement/check_AoP_revision.py` | `aop_supplement/reference_results/` contains the refined tables, `Jost_IR_diagnostic.csv` and `audit_summary.json` |
| Additional scalar stress algebra, 20 independent on-shell two-jets | `aop_supplement/check_scalar_Ward.py` | `aop_supplement/reference_results/scalar_Ward_jets.json` |
| Selection of refined printed rows; no new integration | `aop_supplement/export_publication_tables.py` | `aop_supplement/reference_results/table1_swave_AoP.csv`, `table2_AoP.csv`, `publication_table_sources.json` |

**Refined rows.** Table 1 upper block uses r_min=1e-14, rtol=2e-13,
atol=2e-15. Table 2 uses 96 Gauss–Legendre nodes at the same fine tolerances.
Z=2.5086142772483844. The supplement retains both coarser s-wave runs and the
48-node packet for comparison. The old tables/outputs are not overwritten.
The phase proposition is Proposition 9; Proposition 2 labels the velocity map.

**Comparison rules.** The independent audit uses an absolute tolerance of
`1e-13` only for printed B1/B3 normalization residuals. The symbolic output
and the other numerical fields are compared exactly after elapsed-time suffixes
are removed. For `Jost_IR_diagnostic.csv`, frequencies are compared exactly,
the scaled amplitude uses `rtol=5e-9` and `atol=5e-11`, and both recorded
norm defects must lie in `[0, 1e-6]`. See the
[independent audit](paper_II/independent_audit/README.md) and
[supplementary calculations](paper_II/aop_supplement/README.md).

## Notes
* `num_checks.py` integrates the interior radial equation exactly from the horizon-Jost condition (state variable ln(1-r) near the horizon) and, for K >= 1e3, starts from first-order Liouville-Green data at K eta = 200 (`endpoint_coeffs_LG`), validated against the full integration at K = 1e3.
* In `generic_ray_check.py` only |beta_H| is compared with Eq. (6.6): the overall phase relative to the Jost convention contains exp(i kappa c(lhat)) from the constant offset between the s-wave tortoise time and the CMTV conformal time at the horizon.
* Uncertainties quoted in the papers: C_dyn +/- 4e-10 (dominated by the cutoff test of the fitted tail), Delta_dyn +/- 2e-10 (conservative; the convergence study gives 7e-11), see `convergence_results.txt` and `kasner_complete_rset_output.txt`.


## paper_III — fermionic Weyl/Dirac study

All paths below are relative to `reproducibility/paper_III/`. Checks are independent numerical/symbolic diagnostics and do not replace the proofs in the manuscript. Full mode may be CPU-intensive.

| Result | Source code | Frozen outputs / validation |
|---|---|---|
| Future-time Dirac evolution and modal CAR algebra | `work1/check_time_orientation.py` | `work1/Work1_checks.json`, `work1/Work1_checks.txt` |
| All-order symbol recurrence: elementary Clifford identities and exponent bookkeeping | `work2/checks/check_symbol_algebra.py` | `work2/checks/check_symbol_algebra_output.txt` |
| Invariant spin-harmonic kernel, Jacobi/Bessel scaling, Weyl and local stress diagnostics | `work3/checks/` scripts | `work3/checks/invariant_harmonics_*.{csv,json}`, `continuum_*.{csv,json}` |
| Continuum CAR / conservation / numerical mass cases | `work4/checks/check_work4.py` | `work4/checks/work4_results.json`, `ward_identity_refinement.csv`, `terminal_bounds.csv` |
| Universal transition F, nonmonotone window, high-frequency refinement, C_F | `work5/checks/{F_refinements,convergence_CF,CF_tail_refinement}.py` | `work5/checks/F_refinements.csv`, `CF_quadrature_convergence.csv`, `CF_cutoff_convergence.csv` |
| Fixed-window Hilbert--Schmidt coefficient and toy multiplicity audit | `work6/checks/{window_CF,check_hs_counting}.py` | `work6/checks/window_integrals.csv`, `counting_validation.json` |
| Massive operator transition block vs averaged probability | `work7/checks_work7.py` | `work7/massive_operator_diagnostics.json` |
| Pure CAR state / mixed restriction / rank-one particle-antiparticle algebra | `work8/checks/check_priority4.py` | `work8/results/priority4_exact_CAR_checks.json` |

Run `bash run_quick_checks_paper_III.sh` or `bash run_full_reproduction_paper_III.sh`
from the root. The full runner checks each terminal-bound row, reconstructs
`additive_error` and `bound_ratio`, and records the floating-point roundoff
budgets for those columns. See the [Paper III guide](paper_III/README.md).

## Historical Paper III diagnostics

`paper_III/legacy_phase2/` contains the supplied phase-0 evolution, phase-2 counts,
finite-mode packet and local stress-identity scripts, with their original source
bytes and documented auxiliary sign convention. The Work1–Work8 sources provide
the production checks.
