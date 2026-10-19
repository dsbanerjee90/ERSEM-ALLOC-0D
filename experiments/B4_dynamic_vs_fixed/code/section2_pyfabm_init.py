# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
#
# Section 2 reusable PyFABM initialisation for the B1--B4 model experiments.
# Extracted from the verified ersem_alloc_fixed_fc_sensitivity.py setup.
# ENV and build_model() are preserved, with only the YAML location changed
# to this experiment package's own copied input file.
# This module does not run a model experiment or write to disk.
from pathlib import Path
import numpy as np
import pyfabm

HERE = Path(__file__).resolve().parent
YAML = HERE.parent / 'inputs' / 'fabm.yaml'
FC_REF = 0.24
TEMPERATURE_K = 291.0
ENV = {
    'cell_thickness': 1.0,
    'temperature': TEMPERATURE_K - 273.15,
    'practical_salinity': 35.0,
    'density': 1025.0,
    'pressure': 1.0,
    'N3_flux/flux': 0.0,
    'N4_flux/flux': 0.0,
    'longitude': -4.0,
    'latitude': 50.0,
    'surface_downwelling_shortwave_flux': 50.0,
    'gelbstoff_absorption_satellite': 0.0,
    'wind_speed': 4.0,
    'mole_fraction_of_carbon_dioxide_in_air': 420.0,
    'bottom_stress': 0.0,
    'number_of_days_since_start_of_the_year': 180.0,
}
ENV.update({f'P{i}/fC_state': FC_REF for i in range(1, 5)})


def build_model():
    if not YAML.is_file():
        raise FileNotFoundError(f"Missing copied FABM configuration: {YAML}. Run bash scripts/prepare_scylla.sh first.")
    m = pyfabm.Model(str(YAML))
    unknown = sorted({d.name for d in m.dependencies} - ENV.keys())
    if unknown:
        raise RuntimeError(f'Missing dependency settings: {unknown}')
    for dep in m.dependencies:
        dep.value = ENV[dep.name]
    m.setCellThickness(1.0)
    if not m.checkReady():
        raise RuntimeError('PyFABM reports that model dependencies are not ready')
    names = [s.name for s in m.state_variables]
    ic, ih = names.index('P1/c'), names.index('P1/Chl')
    base = np.array(m.state, dtype=float, copy=True)
    if not np.isfinite(base).all() or base[ic] <= 0:
        raise RuntimeError('Invalid original FABM state / P1 carbon')
    fc = m.findDependency('P1/fC_state')
    sw = m.findDependency('surface_downwelling_shortwave_flux')
    par = m.findDiagnosticVariable('light/parEIR')
    if any(variable is None for variable in (fc, sw, par)):
        raise RuntimeError('Missing P1/fC_state, shortwave or light/parEIR variable')
    return m, base, ic, ih, fc, sw, par
