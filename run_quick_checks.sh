#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT/reproducibility"

echo "[1/3] Symbolic RSET audit"
python check_kasner_rset_symbolics.py

echo "[2/3] Analytic UV-tail audit"
python check_kasner_uv_tails.py

echo "[3/3] Endpoint-Hankel channel audit"
python reproduce_channel_C.py

echo "Quick checks completed."
