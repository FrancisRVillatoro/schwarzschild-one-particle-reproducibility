# Independent audit of Paper II (Annals of Physics version)

Paper II: *Modewise freezing, terminal particle structures and local stress in the Schwarzschild interior* (F. R. Villatoro, submitted to Annals of Physics, 2026).

These two scripts re-derive the symbolic claims of the manuscript and reproduce its Table 1 with an integrator written separately from the production scripts of `reproducibility/paper_II/` (`num_checks.py`, `run_part1.py`, `run_part2.py`, `generic_ray_check.py`, `paper2_symbolic_checks.py`). They do not replace those scripts; they are an independent cross-check of the same results. They were written with the assistance of generative AI tools (ChatGPT, OpenAI; Claude, Anthropic) during the verification of the manuscript, and were reviewed and run by the author, who takes responsibility for them.

## Files

| File | Content | Run time (single core) |
|---|---|---|
| `audit_sym.py` | SymPy re-derivation of the geometry, the pump of Eq. (2.3), the Cortez-type variables of Eq. (3.5), the mass function and the deviations of Lemma 5, the s-wave series (6.2), the radial geodesic expansions (8.9), the Frobenius coefficient of Appendix F, the prefactors of the three particle counts, Eqs. (3.2), (3.4) and (3.6), and the curvature-coupling term of Remark 18 | about 10 s |
| `audit_sym_output.txt` | Output of `audit_sym.py` | |
| `audit_num.py` | (A) Weber--Schafheitlin integral (6.5), the constants $D$ and $A_s$ of Proposition 9, and the Abel value $2/\pi^2$ (mpmath). (B1) Massless s wave, $\kappa|\beta_H|^2$ for $\kappa = 2$ to $128$ and the phase at $\kappa = 512$ (Table 1, top). (B2) Validation of the Liouville--Green start at $K\eta = 200$ against the start at the horizon. (B3) Fixed angular towers $l = 1, 2, 3$, $K = 10^4$ to $10^7$, compared with Table 1 (bottom) and with Theorem 10, Eqs. (6.8)--(6.10). (B4) Directional Born benchmark of Eq. (6.6) along the rays $\hat l^2 = 1/2$ and $4/5$ | about 4 min (about 2.5 min with `SKIP_A=1`) |
| `audit_num_output.txt` | Output of `audit_num.py` | |

## Method of `audit_num.py`

Units $2M = 1$, $\psi = r\phi$, $\eta = -r - \ln(1-r)$ increasing towards the horizon, Jost phase $\psi \sim e^{-iK\eta + i\pi/4}/\sqrt{2K}$. The radial map $r(\eta) = 1 + W_0(-e^{-1-\eta})$ uses the principal branch of the Lambert function. The equation is integrated in $\eta$ from the horizon side ($r_0 = 1 - 10^{-9}$), or from first-order Liouville--Green data at $K\eta = 200$ when $K \ge 10^3$, down to $r = 0.3$, and then in $x = \ln r$ down to $r = 10^{-12}$ (DOP853, relative tolerance $10^{-12}$, absolute tolerance $10^{-16}$). The Frobenius coefficients $(a, b)$ of $\phi = a + b\ln r + \dots$ at $r = 10^{-12}$ are projected onto the endpoint-Hankel pair of the Cortez variables, which gives $(\alpha_H, \beta_H)$; the identity $|\alpha_H|^2 - |\beta_H|^2 = 1$ is monitored and holds to $10^{-11}$ or better in all cases. The printed values labelled "paper" are the entries of Table 1 of the manuscript; those labelled "formula" are Eqs. (6.8)--(6.10).

## How to run

```
python3 audit_sym.py > audit_sym_output.txt
python3 audit_num.py > audit_num_output.txt        # SKIP_A=1 python3 audit_num.py skips part (A)
```

Tested with Python 3.13, NumPy 2.5, SciPy 1.18, SymPy 1.14 and mpmath 1.3. Table 1 is reproduced to all printed digits; the relative deviations of the directional benchmark are approximately 0.5% at $\kappa = 10$ (one value is 0.5013%) and 0.11--0.37% at $\kappa = 40$.

## Release integration notice (v1.2.0)

The submitted output files have been preserved byte-for-byte. The title above has been harmonized with the Annals of Physics TEX in the available working archive, which has a different title from the independent-audit ZIP's addendum. Please confirm the definitive submission title before publishing. The fast-check helper is separate from the supplied numerical audit; it exercises a representative subset without the computationally expensive special-function integrals and horizon-start K=10^4 cross-check. Full verification of those steps must not be inferred from the fast check.
