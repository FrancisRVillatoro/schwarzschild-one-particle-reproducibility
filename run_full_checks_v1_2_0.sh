#!/usr/bin/env bash
set -euo pipefail
unset PYTHONOPTIMIZE
ROOT="$(cd "$(dirname "$0")" && pwd)"
PY="${SCHWARZSCHILD_PYTHON:-python3}"
"$PY" -c 'import sys; sys.exit(0 if sys.version_info >= (3,11) else "Python 3.11 or newer is required")'
exec "$PY" "$ROOT/tools/validate_v120.py" --mode full "$@"
