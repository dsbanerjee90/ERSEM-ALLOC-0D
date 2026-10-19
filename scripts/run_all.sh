#!/usr/bin/env bash
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"
# The same external PyFABM path used successfully in the original 0D work.
EXTERNAL_PYFABM="${PYFABM_PYTHONPATH:-/work/dba/ERSEM-ALLOC-0D/paper/pyfabm_tests/python}"
if [ ! -d "$EXTERNAL_PYFABM/pyfabm" ]; then
  echo "STOP: external PyFABM not found: $EXTERNAL_PYFABM" >&2; exit 2
fi
export PYTHONPATH="$EXTERNAL_PYFABM${PYTHONPATH:+:$PYTHONPATH}"
for file in \
  results/geider1997_ersem_balanced_chlc_model.csv \
  results/ersem_alloc_fixed_fc_sensitivity.csv \
  results/ersem_alloc_p1_photosynthesis_irradiance.csv; do
  if [ -e "$file" ]; then
    echo "STOP: a previous result exists: $file" >&2
    echo 'Keep it for provenance. For a new run, work in a fresh checkout or manually archive generated results.' >&2
    exit 3
  fi
done
python scripts/bootstrap_scylla.py
python scripts/smoke_test.py > results/smoke_test.log 2>&1 || { tail -50 results/smoke_test.log; exit 4; }
echo 'PASS PyFABM smoke test (full log: results/smoke_test.log)'
echo '[1/3] Geider reference balance'
python scripts/run_geider_reference.py > results/run_geider_reference.log 2>&1 || { tail -80 results/run_geider_reference.log; exit 5; }
echo '[2/3] Prescribed fC Chl:C sensitivity'
python scripts/run_balanced_fc_sweep.py > results/run_balanced_fc_sweep.log 2>&1 || { tail -80 results/run_balanced_fc_sweep.log; exit 6; }
echo '[3/3] P1 gross fixation P-I sensitivity'
python scripts/run_pi_sweep.py > results/run_pi_sweep.log 2>&1 || { tail -80 results/run_pi_sweep.log; exit 7; }
python scripts/verify_against_reference.py
python scripts/plot_supp_geider_pi.py
printf '\nDONE: files in results/ and figures/.\n'
ls -lh figures/supp_geider_pi.{png,pdf,svg}
