# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
# Reproduce a Figure 12-style ERSEM balanced-growth Chl:C comparison using
# published nutrient-replete Thalassiosira pseudonana culture observations.
# Observations: Geider et al. (1997), Figure 3A; original data Geider (1984).
# Methodological precedent: Butenschoen et al. (2016), Figure 12.
# Existing ERSEM/FABM code, simulations, and YAML files are never modified.

from pathlib import Path
import csv

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyfabm

ROOT = Path(__file__).resolve().parents[1]  # This new PAPER workspace
YAML = ROOT / 'configs/fabm.yaml'
OUT = ROOT / 'results'
OBS = ROOT / 'observations/geider1997_thalassiosira_observations.csv'
OUT.mkdir(parents=True, exist_ok=True)

# Geider et al. (1997), Table 2: 291 K, continuous illumination, nutrient-replete.
TEMPERATURE_C = 291.0 - 273.15
FC_REFERENCE = 0.24  # Recovers reference carbon-fixation machinery capacity.

# This approximate spectral conversion requires an assumed PAR photon spectrum.
# 4.57 umol photons/J is a conventional approximation for broadband PAR.
# For exact spectral conversion, measured/incidental spectral irradiance is needed.
PAR_UMOL_PER_J = 4.57

# Illuminated continuously at prescribed intensities, not a time-varying day/night cycle.
# The model is generic P1 diatoms, not a species-specific T. pseudonana culture fit.
SHORTWAVE_GRID = np.linspace(2.0, 350.0, 48)  # W m-2, broadband forcing
CHLC_GRID = np.linspace(0.002, 0.10, 221)     # g Chl a per g C (equivalent mg/mg)

# Keep Figure 12 comparison simple: one model curve and the published observations.
# Set to True only to add an explanatory, non-fitted machinery-allocation sensitivity.
INCLUDE_FC_SENSITIVITY = False
FC_CASES = [FC_REFERENCE] + ([0.38] if INCLUDE_FC_SENSITIVITY else [])

if not OBS.is_file():
    raise FileNotFoundError(f'Missing digitised observations: {OBS}')
with OBS.open(newline='') as f:
    rr = list(csv.DictReader(f))
if len(rr) != 11:
    raise ValueError(f'Expected 11 digitised Fig. 3A observations, found {len(rr)}')
x_obs = np.asarray([float(r['PAR_umol_photons_m2_s']) for r in rr])
y_obs = np.asarray([float(r['ChlC_g_g']) for r in rr])
if not (np.isfinite(x_obs).all() and np.isfinite(y_obs).all()
        and (x_obs >= 0).all() and (y_obs > 0).all()):
    raise ValueError('Nonphysical or nonfinite digitised observations')

m = pyfabm.Model(str(YAML))
env = {
    'cell_thickness': 1.0,
    'temperature': TEMPERATURE_C,
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
for pft in range(1, 5):
    env[f'P{pft}/fC_state'] = FC_REFERENCE

for dep in m.dependencies:
    dep.value = env[dep.name]
m.setCellThickness(1.0)
assert m.checkReady(), 'FABM model not ready'

names = [v.name for v in m.state_variables]
ic = names.index('P1/c')
ichl = names.index('P1/Chl')
initial = np.array(m.state, dtype=float, copy=True)
c0 = float(initial[ic])
if c0 <= 0:
    raise ValueError('Initial P1 carbon must be positive')
fc_dep = m.findDependency('P1/fC_state')
sw_dep = m.findDependency('surface_downwelling_shortwave_flux')
par = m.findDiagnosticVariable('light/parEIR')  # W m-2 of PAR
if par is None:
    raise RuntimeError('Missing light/parEIR: cannot compare irradiance units')


def balanced_residual(theta, f_c, sw, return_par=False):
    """d ln(Chl/C) / dt in per-day units, keeping other state variables fixed."""
    m.state[:] = initial
    m.state[ichl] = theta * c0
    fc_dep.value = float(f_c)
    sw_dep.value = float(sw)
    rates = np.asarray(m.getRates(), dtype=float)
    if not np.isfinite(rates).all():
        raise RuntimeError(f'Nonfinite rates at fC={f_c}, SW={sw}, theta={theta}')
    residual = 86400.0 * (rates[ichl] / (theta*c0) - rates[ic] / c0)
    if return_par:
        return float(residual), float(par.value)
    return float(residual)


def solve_chlc(f_c, sw):
    vals = [balanced_residual(t, f_c, sw) for t in CHLC_GRID]
    roots = []
    for i in range(len(CHLC_GRID)-1):
        a, b = vals[i:i+2]
        if abs(a) < 1e-10:
            roots.append(float(CHLC_GRID[i]))
        elif a*b < 0:
            roots.append(float(brentq(
                lambda t: balanced_residual(t, f_c, sw),
                float(CHLC_GRID[i]), float(CHLC_GRID[i+1]), xtol=1e-12
            )))
    if abs(vals[-1]) < 1e-10:
        roots.append(float(CHLC_GRID[-1]))
    roots = sorted(set(round(v, 11) for v in roots))
    if len(roots) != 1:
        return np.nan, np.nan, np.nan, ('NO_ROOT' if not roots else 'MULTIPLE_ROOTS')
    theta = float(roots[0])
    residual, model_par = balanced_residual(theta, f_c, sw, return_par=True)
    return theta, model_par, residual, 'ONE_ROOT'


results=[]
for f_c in FC_CASES:
    for sw in SHORTWAVE_GRID:
        theta, par_w, residual, status = solve_chlc(f_c, float(sw))
        results.append({
            'fC': f_c,
            'shortwave_W_m2': float(sw),
            'model_PAR_W_m2': par_w,
            'model_PAR_umol_photons_m2_s': par_w*PAR_UMOL_PER_J,
            'balanced_ChlC_g_g': theta,
            'balance_residual_day': residual,
            'status': status,
        })
    print(f'fC={f_c:.2f}: {sum(r["status"]=="ONE_ROOT" for r in results if r["fC"]==f_c)}/{len(SHORTWAVE_GRID)} single-root cases')

out_csv = OUT / 'geider1997_ersem_balanced_chlc_model.csv'
with out_csv.open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(results[0]))
    writer.writeheader()
    writer.writerows(results)

fig, ax = plt.subplots(figsize=(7.5, 5.25), constrained_layout=True)
for f_c in FC_CASES:
    rr = sorted([r for r in results if r['fC'] == f_c and r['status'] == 'ONE_ROOT'],
                key=lambda r: r['model_PAR_umol_photons_m2_s'])
    if rr:
        ax.plot([r['model_PAR_umol_photons_m2_s'] for r in rr],
                [r['balanced_ChlC_g_g'] for r in rr], lw=2.2,
                linestyle=('-' if f_c == FC_REFERENCE else '--'),
                label=(f'ERSEM-ALLOC (reference $f_C={f_c:.2f}$)' if f_c == FC_REFERENCE
                       else f'Fixed $f_C={f_c:.2f}$ sensitivity'))
ax.scatter(x_obs, y_obs, color='black', s=38, marker='o', zorder=4,
           label=r'$T.$ pseudonana observations (Geider et al., 1997; digitised)')
ax.set(xlim=(0,600), ylim=(0,0.065),
       xlabel=r'Irradiance ($\mu$mol photons m$^{-2}$ s$^{-1}$)',
       ylabel=r'Chlorophyll $a$:carbon (g Chl $a$ g$^{-1}$ C)')
ax.grid(alpha=0.15)
ax.legend(frameon=False, fontsize=8.8, loc='upper right')
ax.text(0.02, 0.03,
        '291 K, continuous illumination; culture observations digitised from Fig. 3A',
        transform=ax.transAxes, fontsize=8, color='0.32')
fig.savefig(OUT / 'ersem_alloc_geider_figure12.png', dpi=250)
fig.savefig(OUT / 'ersem_alloc_geider_figure12.pdf')
plt.close(fig)

print('Observations: 11 digitised points, Geider et al. (1997), Fig. 3A')
print('Model: P1 generic diatoms, fixed initial N/P, fC=0.24, T=291 K')
print('PAR conversion: model W m-2 x 4.57 umol photons J-1 (approximate)')
print('Not a direct experimental validation of the new fC optimiser')
print('Wrote:', out_csv)
print('Wrote:', OUT / 'ersem_alloc_geider_figure12.png')
print('Wrote:', OUT / 'ersem_alloc_geider_figure12.pdf')
