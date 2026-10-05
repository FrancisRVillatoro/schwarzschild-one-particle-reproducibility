# Reproducibility package for the two Schwarzschild-interior papers

Prepared 4 October 2026. Python 3 with numpy, scipy, sympy, mpmath. Units: the papers' conventions
(paper I: G = c = hbar = 1, Kasner coordinates with y = 1; paper II: 2M = 1 unless stated).
Each script writes to standard output; the file `<script>_output.txt` next to it is the output
of the run used for the manuscripts. Single-core run times are indicated.

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

## paper_II  (Particle number, terminal states and local stress near the Schwarzschild singularity)

| Manuscript item | Script | Output file | Run time |
|---|---|---|---|
| Table 1 (top): massless s wave, kappa |beta_H|^2 for kappa = 2 ... 128 | `num_checks.py` | `num_checks_output.txt` | 1 min |
| Proposition 6 phase check, beta_H sqrt(kappa) -> D(1+i); Table 1 (bottom), fixed angular towers l = 1,2,3, K = 1e2 ... 1e7, against Theorem 7 | `run_part1.py` (needs `num_checks.py`) | `run_part1_output.txt` | 4 min |
| Table 2: r^3 Delta rho_E /((3/2)|B_*|^2) for a conditioned packet (l = 0; l = 1 pole; l = 1 equator) and E = 0, 0.5, 2 | `run_part2.py` (needs `num_checks.py`) | `run_part2_output.txt` | 3 min |
| Eq. (6.6), generic rays lhat^2 = 1/2 and 4/5: |beta_H|^2 kappa^4 against 1/(36 pi^2 lhat^4) | `generic_ray_check.py` (needs `num_checks.py`) | `generic_ray_check_output.txt` | 1 min |
| Eq. (2.4) for Q, H and its zero at r = 3/4 (Proposition 3), velocity map normalization, sign of the Hankel normalization (2 Im(Pbar Q) = -pi/2), Weber-Schafheitlin integral, D and A_s, Abel-regularized integral 2/pi^2, Theorem 7 algebra, excitation coupling term of Remark 13 | `paper2_symbolic_checks.py` | `paper2_symbolic_checks_output.txt` | 3 min |
| Summary of the numerical checks as used in the text | | `numerical_checks_results.txt` | |

## Notes
* `num_checks.py` integrates the interior radial equation exactly from the horizon-Jost condition (state variable ln(1-r) near the horizon) and, for K >= 1e3, starts from first-order Liouville-Green data at K eta = 200 (`endpoint_coeffs_LG`), validated against the full integration at K = 1e3.
* In `generic_ray_check.py` only |beta_H| is compared with Eq. (6.6): the overall phase relative to the Jost convention contains exp(i kappa c(lhat)) from the constant offset between the s-wave tortoise time and the CMTV conformal time at the horizon.
* Uncertainties quoted in the papers: C_dyn +/- 4e-10 (dominated by the cutoff test of the fitted tail), Delta_dyn +/- 2e-10 (conservative; the convergence study gives 7e-11), see `convergence_results.txt` and `kasner_complete_rset_output.txt`.
