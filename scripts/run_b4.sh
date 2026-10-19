#!/usr/bin/env bash
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
# Optional: reproduce main manuscript Fig. 3 B4 0D dynamics from this PAPER directory.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"
export PYTHONPATH="${PYFABM_PYTHONPATH:-/work/dba/ERSEM-ALLOC-0D/paper/pyfabm_tests/python}${PYTHONPATH:+:$PYTHONPATH}"
python scripts/bootstrap_scylla.py
EXP="$HERE/experiments/B4_dynamic_vs_fixed"
OUT="${B4_OUTPUT_DIR:-$EXP/results/B4_dt150_fresh}"
LOG="$HERE/results/$(basename "$OUT").log"
if [ -e "$OUT" ]; then echo "STOP: B4 fresh output exists, no overwrite: $OUT" >&2; exit 3; fi
mkdir -p "$(dirname "$LOG")"
(
 cd "$EXP"
 python code/B4_dynamic_vs_fixed.py --days 7 --dt 150 --lights 20,100 --output "$OUT"
) > "$LOG" 2>&1 || { tail -80 "$LOG"; exit 4; }
python "$EXP/figures/plot_fig3.py" --input "$OUT" --output "$EXP/figures"
ls -lh "$EXP/figures/fig3.pdf" "$EXP/figures/fig3.png"
