# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
# Compare already computed B4 summary metrics at 150-s and 75-s steps.
from __future__ import annotations
import argparse, csv
from pathlib import Path

def load(path):
    with Path(path).open(newline='') as f: rows=list(csv.DictReader(f))
    return {(float(r['shortwave_W_m2']),r['treatment']):r for r in rows}

p=argparse.ArgumentParser()
p.add_argument('--coarse',type=Path,required=True)
p.add_argument('--fine',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
coarse=load(a.coarse/'B4_paired_comparison_summary.csv')
fine=load(a.fine/'B4_paired_comparison_summary.csv')
keys=('fC_dynamic_final', 'P1_C_percent_change_dynamic_vs_fixed',
      'P1_gross_C_fixation_specific_pct_change','total_inorganic_N_specific_pct_change',
      'P1_PO4_uptake_specific_pct_change')
rows=[]
for k in sorted(coarse):
    if k not in fine: raise RuntimeError('Missing fine timestep case: '+str(k))
    for metric in keys:
        x,y=float(coarse[k][metric]),float(fine[k][metric])
        rows.append({'shortwave_W_m2':k[0],'treatment':k[1],'metric':metric,
                     'dt150':x,'dt75':y,'abs_difference':abs(x-y)})
if a.output.exists(): raise FileExistsError(a.output)
a.output.parent.mkdir(parents=True,exist_ok=True)
with a.output.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
for k in sorted(coarse):
    relevant=[r for r in rows if (r['shortwave_W_m2'],r['treatment'])==k]
    print(k, 'max metric difference=',max(r['abs_difference'] for r in relevant),
          '(note: metrics have mixed units)')
print('Saved',a.output)
