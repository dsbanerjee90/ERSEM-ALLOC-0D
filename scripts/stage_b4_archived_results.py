#!/usr/bin/env python3
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
"""Copy original validated B4 CSV baseline into selected PAPER workspace without changes."""
from pathlib import Path
import shutil,os
ROOT=Path(__file__).resolve().parents[1]
SOURCE=Path(os.environ.get('ERSEM_ALLOC_ROOT','/work/dba/ERSEM-ALLOC-0D'))
ORIG=SOURCE/'paper/ERSEM_ALLOC_Section2_repro_v1/results/B4_dt150'
DEST=ROOT/'experiments/B4_dynamic_vs_fixed/reference/B4_dt150_original'
files=sorted(ORIG.glob('B4_*.csv'))
if not files or not (ORIG/'B4_paired_comparison_summary.csv').exists():
    raise SystemExit('No B4 CSV files found: '+str(ORIG))
if DEST.exists() and any(DEST.iterdir()):
    raise SystemExit('STOP: B4 baseline already staged; will not overwrite: '+str(DEST))
DEST.mkdir(parents=True,exist_ok=True)
for f in files:shutil.copy2(f,DEST/f.name)
print('Staged',len(files),'original B4 CSVs in',DEST)
print('This is an ARCHIVED baseline, not a fresh B4 run.')
