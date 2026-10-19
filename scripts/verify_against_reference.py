#!/usr/bin/env python3
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
"""Numerically compare new experiment CSVs with archived original outputs; no refitting."""
import csv, json, hashlib, math, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
FILES=['geider1997_ersem_balanced_chlc_model.csv',
       'ersem_alloc_fixed_fc_sensitivity.csv',
       'ersem_alloc_p1_photosynthesis_irradiance.csv']
RTOL=2e-6
ATOL=1e-8

def load(p):
    with p.open(newline='') as f:
        rd=csv.DictReader(f);return list(rd.fieldnames or []),list(rd)
def verify():
    report={'relative_tolerance':RTOL,'absolute_tolerance':ATOL,'comparisons':{}}
    for name in FILES:
        reference=ROOT/'reference'/name
        current=ROOT/'results'/name
        if not reference.is_file() or not current.is_file():
            raise RuntimeError('Missing comparison CSV: '+str(current if not current.is_file() else reference))
        refcols, rr=load(reference); cols, rows=load(current)
        if refcols!=cols:raise RuntimeError(f'{name}: columns differ. Original={refcols}, new={cols}')
        if len(rr)!=len(rows):raise RuntimeError(f'{name}: row count differs, {len(rr)} vs {len(rows)}')
        max_abs=0.0;max_rel=0.0;ndiff=0
        for i,(a,b) in enumerate(zip(rr,rows)):
            for k in cols:
                try: av=float(a[k]);bv=float(b[k]); numeric=True
                except (ValueError,TypeError):numeric=False
                if not numeric:
                    if a[k]!=b[k]:raise RuntimeError(f'{name} row {i} {k}: {a[k]!r} != {b[k]!r}')
                    continue
                if math.isnan(av) and math.isnan(bv):continue
                if not math.isfinite(av) or not math.isfinite(bv):
                    if av!=bv:raise RuntimeError(f'{name} row {i} {k}: nonfinite mismatch')
                    continue
                err=abs(av-bv)
                max_abs=max(max_abs,err)
                max_rel=max(max_rel,err/max(abs(av),abs(bv),1e-30))
                if err>ATOL+RTOL*max(abs(av),abs(bv)):
                    ndiff+=1
                    if ndiff<3: print(f'MISMATCH: {name} row={i} {k} reference={av:.10g} new={bv:.10g}')
        report['comparisons'][name]={ 'rows':len(rows),'max_abs_difference':max_abs,
                                      'max_relative_difference':max_rel,'mismatch_values':ndiff,
                                      'passed':ndiff==0}
        print(f'{"PASS" if ndiff==0 else "FAIL"}: {name}, rows={len(rows)}, max_abs={max_abs:.4g}, mismatches={ndiff}',flush=True)
    p=ROOT/'results/reference_comparison.json';p.write_text(json.dumps(report,indent=2)+'\n')
    if any(not v['passed'] for v in report['comparisons'].values()):
        raise SystemExit('FAIL: baseline numerical disagreement. Do not interpret the new figure as reproduced. Inspect PyFABM binary, YAML and model version.')
    print('PASS: all regenerated CSVs agree numerically with archived baselines.')
if __name__=='__main__':verify()
