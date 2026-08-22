#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "$SCRIPT_DIR/run_pipeline.py" \
  --out "$SCRIPT_DIR/data/practical_50_loci" \
  --taxa 51 \
  --loci 50 \
  --length 500 \
  --seed 6406001 \
  "$@"

python3 "$SCRIPT_DIR/evaluate.py" \
  --run "$SCRIPT_DIR/data/practical_50_loci"
