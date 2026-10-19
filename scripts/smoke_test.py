#!/usr/bin/env python3
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
"""Read-only PyFABM initialization/diagnostic availability smoke test."""
from pathlib import Path
import sys
root=Path(__file__).resolve().parents[1]
config=root/'configs/fabm.yaml'
if not config.is_file():
    raise SystemExit('Missing configs/fabm.yaml: first run python scripts/bootstrap_scylla.py')
import numpy as np
import pyfabm
print('Python:',sys.executable,flush=True)
print('PyFABM:',pyfabm.__file__,flush=True)
m=pyfabm.Model(str(config))
for name in ('P1/fC_state','surface_downwelling_shortwave_flux'):
    dep=m.findDependency(name)
    if dep is None: raise RuntimeError('Missing dependency '+name)
for name in ('P1/fO3PIc','light/parEIR'):
    diag=m.findDiagnosticVariable(name)
    if diag is None:raise RuntimeError('Missing diagnostic '+name)
for name in ('P1/c','P1/Chl'):
    if name not in [v.name for v in m.state_variables]:
        raise RuntimeError('Missing state '+name)
print('PASS: P1 states, fC dependency, irradiance and gross-fixation diagnostics present.',flush=True)
