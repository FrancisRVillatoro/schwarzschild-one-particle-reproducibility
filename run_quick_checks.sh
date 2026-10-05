#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"

echo "[1/3] Paper I: adiabatic symbolic audit"
(cd "$ROOT/reproducibility/paper_I" && python kasner_adiabatic_symbolics.py)

echo "[2/3] Paper I: equal-time point splitting"
(cd "$ROOT/reproducibility/paper_I" && python kasner_point_splitting.py)

echo "[3/3] Paper II: symbolic audit"
(cd "$ROOT/reproducibility/paper_II" && python paper2_symbolic_checks.py)

echo "Quick checks completed."
