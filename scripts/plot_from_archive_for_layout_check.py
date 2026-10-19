#!/usr/bin/env python3
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
"""Only for layout inspection without PyFABM. Never claim this is a fresh run.

Run on a new clone with --reference-only to COPY archived CSVs into results.
It refuses to overwrite generated CSVs. For scientific regeneration run run_all.sh.
"""
from pathlib import Path
import shutil,sys
ROOT=Path(__file__).resolve().parents[1]
if sys.argv[1:]!=['--reference-only']:
    raise SystemExit('Only use: python scripts/plot_from_archive_for_layout_check.py --reference-only')
for name in ('ersem_alloc_fixed_fc_sensitivity.csv','ersem_alloc_p1_photosynthesis_irradiance.csv'):
    dest=ROOT/'results'/name
    if dest.exists():raise SystemExit(f'STOP: {dest} exists; will not overwrite.')
    shutil.copy2(ROOT/'reference'/name,dest)
print('REFERENCE-ONLY layout check: NO NEW ERSEM MODEL INTEGRATION OR RATES COMPUTED.')
from plot_supp_geider_pi import main
main()
