# Technical reproduction notes

This directory contains the numerical and symbolic checks supporting the results associated with *Fate of a one-particle state approaching the Schwarzschild spacelike singularity* and the overlapping Kasner fixed-point calculations used in *Universal Quantum Kasner Fixed Point at the Schwarzschild Singularity*.

No manuscript source or PDF is included in this repository.

## 1. Kasner vacuum polarization

Run

```bash
python reproducibility/reproduce_kasner_phi2.py
```

The calculation evolves the terminal Kasner modes independently of the Schwarzschild mode campaign and evaluates the second-order adiabatically subtracted coincidence integral. The expected result is

- `C_dyn ~= -2.62616e-4`,
- `C_IR = 1/(240*pi^2)`,
- `C_K ~= 1.59555e-4`,
- `q_K = (9/2) C_K ~= 7.1800e-4`.

`check_kasner_massive_IR.py` checks explicitly that the shrinking exact massive-mode layer leaves no additional finite residue.

## 2. Complete terminal Kasner RSET

The production calculation is

```bash
python reproducibility/kasner_complete_rset.py
```

with numerical routines in `kasner_rset_core.py`. The production path uses the regulator-free massless anisotropy integral `integrate_delta_massless_interval`; finite-mass full-component routines retained in the core file are diagnostics and are not used to quote the final tensor.

Expected central values are

- `Delta_dyn ~= 7.0609e-5` (quoted uncertainty `3e-9`),
- `Delta_IR = -523/(498960*pi^2) ~= -1.0620286004e-4`,
- `Delta_K ~= -3.5594e-5`,
- `C_rho ~= 5.9211e-5`,
- `C_parallel ~= 1.3018e-4`,
- `C_perp ~= 1.6577e-4`,
- Schwarzschild coefficients `(lambda_t,lambda_r,lambda_perp) ~= (2.6361,-1.1990,3.3568)e-3`.

The symbolic audit

```bash
python reproducibility/check_kasner_rset_symbolics.py
```

rederives the fourth-order angular subtraction, verifies `J_0(lambda)=-1/(63 lambda^3)`, and evaluates the nonuniform residue exactly as `-523/(498960*pi^2)`. The sixth-order ultraviolet coefficients used to control the tails are

- `c7 = -1/2835` for the vacuum-polarization angular bracket,
- `d7 = 4/2835` for the anisotropic RSET angular bracket.

Run

```bash
python reproducibility/check_kasner_uv_tails.py
```

to evaluate their closed-form leading tail contributions. The production RSET script uses the exact `d7` coefficient to control the leading ultraviolet tail, estimates the subleading `lambda^-9` term from the stable panel, and assigns a conservative absolute uncertainty to `Delta_dyn`.

## 3. Endpoint-Hankel channel audit

Run

```bash
python reproducibility/reproduce_channel_C.py
```

This reconstructs the direct endpoint-Hankel projections of the four historical channels and the channel-C near-cancellation. The supplied CSV/TXT files record endpoint convergence, preparation sensitivity, and the isolated zero-squeeze matching resonance.

## 4. Interpretation

The numerical values `C_K` and `Delta_K` are reproducible numerical evaluations of convergent renormalized Kasner integrals. The associated analytic statements concern the local anisotropic Schwarzschild-to-Kasner blow-up, the common Unruh/Hartle--Hawking terminal kernel, the fourth-order subtraction structure, conservation/trace reconstruction, and the absence of an additional finite massless-layer residue. Jiang--Jiang values are used only as external benchmarks, not as inputs.
