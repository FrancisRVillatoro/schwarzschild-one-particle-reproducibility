# Paper II — supplementary scalar calculations

Paper: *Particle creation, terminal states and local stress of a scalar field in the Schwarzschild interior*, Francisco R. Villatoro (2026).
These additions do not replace or edit the five production scripts of v1.1.0.
Equation and table labels follow the Paper II manuscript.

## Exact source names and outputs

- `check_AoP_revision.py` performs
  29 algebraic/numerical checks and regenerates `swave_refinement.csv`,
  `normalized_packet.csv`, `packet_refinement.csv`, `Jost_IR_diagnostic.csv`,
  `squeezed_pair.csv`, and `audit_summary.json` in `results/`.
- `check_scalar_Ward.py`: 20 independently constructed on-shell scalar two-jet
  checks; creates `results/scalar_Ward_jets.json`. This is not a duplicate of
  `paper2_symbolic_checks.py` and does not calculate exterior scattering.
- `export_publication_tables.py`: selects and exports the refined rows as
  `table1_swave_AoP.csv`, `table2_AoP.csv`, and `publication_table_sources.json`.
  It computes no new field dynamics and does not refit the original results.

The supplementary scalar script imports `num_checks.py` from the parent
`paper_II/` directory, using the v1.1.0 numerical implementation.

## Traceability of the printed tables

Table 1, upper block: select r_min=1e-14, rtol=2e-13, atol=2e-15 from
`swave_refinement.csv`. Coarser runs are r_min=1e-10, rtol=1e-11, atol=1e-13,
and r_min=1e-12, rtol=2e-12, atol=2e-14. The fixed-tower lower block continues
using the unchanged `../run_part1.py`; the independently written
`../independent_audit/audit_num.py` verifies its displayed numbers.

Table 2: select N_omega=96, rtol=2e-13, atol=2e-15 in `normalized_packet.csv`.
The other runs use 48 nodes with (1e-11,1e-13) and 96 nodes with
(2e-12,2e-14). The horizon-flux normalization is
Z=2.5086142772483844. The nonzero-coefficient ratios are normalization-invariant;
the absolute angular-node column must be divided by Z relative to the raw
v1.1.0 diagnostic. This is NOT a numerical exterior Unruh-transmission calculation.
The legacy defaults are not presented as producing all revised printed digits.

## Run without modifying frozen reference outputs

From the repository root:

```sh
python tools/validate_v120.py --paper II --mode full --report /tmp/v120_paperII.json
```

This runs in a temporary copy and compares newly generated results with the
recorded references. `reference_results/` is not overwritten. The five old
production files are imported, not copied into this new directory.
For a direct local diagnostic, run `check_AoP_revision.py` first (it creates
`results/`), then `check_scalar_Ward.py`, then `export_publication_tables.py`.

The joint report records the environment, regenerated table comparisons and
Ward checks. In `Jost_IR_diagnostic.csv`, it compares frequencies exactly and
the scaled amplitude with `rtol=5e-9` and `atol=5e-11`. It checks both
`norm_error` values against the interval `[0, 1e-6]`, accounting for the
cancellation in the unit Wronskian.
