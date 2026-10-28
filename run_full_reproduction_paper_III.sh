#!/usr/bin/env bash
ROOT="$(cd "$(dirname "$0")" && pwd)"
python "$ROOT/reproducibility/paper_III/run_checks.py" --mode full "$@"
