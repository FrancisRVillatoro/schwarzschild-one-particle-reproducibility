#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT/reproducibility"

echo "[1/6] Kasner vacuum polarization"
python reproduce_kasner_phi2.py

echo "[2/6] Massive-IR consistency check"
python check_kasner_massive_IR.py

echo "[3/6] Complete terminal Kasner RSET"
python kasner_complete_rset.py

echo "[4/6] Symbolic RSET audit"
python check_kasner_rset_symbolics.py

echo "[5/6] Analytic UV-tail audit"
python check_kasner_uv_tails.py

echo "[6/6] Endpoint-Hankel channel audit"
python reproduce_channel_C.py

echo "Full reproduction completed."
