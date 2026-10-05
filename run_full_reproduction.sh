#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"

echo "[1/15] Paper I: vacuum polarization baseline"
(
    cd "$ROOT/reproducibility/paper_I"
    python reproduce_kasner_phi2.py
)

echo "[2/15] Paper I: C_dyn convergence study"
(
    cd "$ROOT/reproducibility/paper_I"
    rm -f conv_C.txt
    python convergence_study.py C \
        baseline target120 Nu32 cut6 n24 rtol3e-11 fit8-14 Nu40 target240
)

echo "[3/15] Paper I: Delta_dyn convergence study"
(
    cd "$ROOT/reproducibility/paper_I"
    rm -f conv_D.txt
    python convergence_study.py D \
        baseline cut8 Nu34 Nl26 target240 rtol1e-12
)

echo "[4/15] Paper I: direct ultraviolet tail fit"
(
    cd "$ROOT/reproducibility/paper_I"
    python tailfit.py
)

echo "[5/15] Paper I: massive-IR check"
(
    cd "$ROOT/reproducibility/paper_I"
    python check_kasner_massive_IR.py
)

echo "[6/15] Paper I: final Kasner RSET"
(
    cd "$ROOT/reproducibility/paper_I"
    python kasner_complete_rset.py
)

echo "[7/15] Paper I: adiabatic symbolic audit"
(
    cd "$ROOT/reproducibility/paper_I"
    python kasner_adiabatic_symbolics.py
)

echo "[8/15] Paper I: equal-time point splitting"
(
    cd "$ROOT/reproducibility/paper_I"
    python kasner_point_splitting.py
)

echo "[9/15] Paper I: legacy RSET symbolic cross-check"
(
    cd "$ROOT/reproducibility/paper_I"
    python check_kasner_rset_symbolics.py
)

echo "[10/15] Paper I: legacy analytic UV-tail cross-check"
(
    cd "$ROOT/reproducibility/paper_I"
    python check_kasner_uv_tails.py
)

echo "[11/15] Paper II: principal numerical checks"
(
    cd "$ROOT/reproducibility/paper_II"
    python num_checks.py
)

echo "[12/15] Paper II: s-wave and fixed-l checks"
(
    cd "$ROOT/reproducibility/paper_II"
    python run_part1.py
)

echo "[13/15] Paper II: generic-ray check"
(
    cd "$ROOT/reproducibility/paper_II"
    python generic_ray_check.py
)

echo "[14/15] Paper II: one-particle stress check"
(
    cd "$ROOT/reproducibility/paper_II"
    python run_part2.py
)

echo "[15/15] Paper II: symbolic audit"
(
    cd "$ROOT/reproducibility/paper_II"
    python paper2_symbolic_checks.py
)

echo "Full reproduction completed."
