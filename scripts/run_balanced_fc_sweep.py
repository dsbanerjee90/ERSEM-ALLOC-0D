# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
#
# ERSEM-ALLOC P1 diatom test harness: sensitivity of balanced-growth
# chlorophyll a:carbon to PRESCRIBED carbon-fixation machinery fraction f_C.
# Methodology follows the PyFABM / Butenschoen et al. (2016) Figure 12 style.
# Observations: T. pseudonana, Geider et al. (1997), Fig. 3A, digitised
# independently in WebPlotDigitizer (original culture data: Geider, 1984).
#
# This is a conditional physiological sensitivity experiment; it does NOT
# solve for optimal f_C or evolve f_C prognostically. Existing YAML, Fortran,
# model runs, previously generated CSVs and figures are NEVER modified.

from pathlib import Path
import csv

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import pyfabm

ROOT = Path(__file__).resolve().parents[1]  # This new PAPER workspace
YAML = ROOT / 'configs/fabm.yaml'
OUTPUT = ROOT / 'results'
OBS = ROOT / 'observations/geider1997_thalassiosira_observations.csv'
VERIFIED_REFERENCE = OUTPUT / 'ersem_reference_vs_legacy_consistency.csv'

FC_REF = 0.24
FC_VALUES = (0.12, 0.24, 0.38, 0.50)
FLH_FIXED = 0.38  # From the allocation formulation; f_V=1-f_LH-f_C.
TEMPERATURE_K = 291.0
PAR_UMOL_PER_J = 4.57   # APPROXIMATE conversion for broadband PAR.
SW_GRID = np.linspace(2.0, 350.0, 48)  # W m^-2: same as verified control.
CHLC_GRID = np.linspace(0.002, 0.100, 221)  # g Chl a / g C.
BALANCE_TOL_PER_DAY = 1e-7

# Model-wide dependencies: all PFTs other than P1 stay at reference f_C.
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


def read_observations():
    if not OBS.is_file():
        raise FileNotFoundError(f'WebPlotDigitizer observations not found: {OBS}')
    with OBS.open(newline='') as handle:
        items = list(csv.DictReader(handle))
    if len(items) != 11:
        raise ValueError(f'Expected 11 Geider Fig. 3A observations, got {len(items)}')
    x = np.array([float(item['PAR_umol_photons_m2_s']) for item in items])
    y = np.array([float(item['ChlC_g_g']) for item in items])
    if (not np.isfinite(x).all() or not np.isfinite(y).all()
            or (x < 0).any() or (y <= 0).any()):
        raise ValueError('Observations contain nonfinite/nonphysical values')
    return x, y


def build_model():
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


def evaluate(ctx, theta, f_c, shortwave):
    """Chl:C relative-growth balance (day^-1) at a controlled instant.

    Only P1 chlorophyll state and fC dependency are varied. All other state
    variables, including nutrient pools, reset to their identical initial
    values at EVERY evaluation. Hence not a full nutrient-quota equilibrium.
    """
    m, base, ic, ih, fc, sw, par = ctx
    m.state[:] = base
    m.state[ih] = float(theta) * float(base[ic])
    fc.value = float(f_c)
    sw.value = float(shortwave)
    rates = np.asarray(m.getRates(), dtype=float)
    if not np.isfinite(rates).all():
        raise RuntimeError(f'Nonfinite rates: fC={f_c}, SW={shortwave}, ChlC={theta}')
    residual = 86400.0 * (rates[ih] / m.state[ih] - rates[ic] / m.state[ic])
    par_wm2 = float(par.value)
    if not np.isfinite(residual) or not np.isfinite(par_wm2) or par_wm2 < 0:
        raise RuntimeError(f'Invalid balance/PAR: fC={f_c}, SW={shortwave}')
    return float(residual), par_wm2


def solve_balanced_chlc(ctx, f_c, shortwave):
    vals = np.asarray([evaluate(ctx, t, f_c, shortwave)[0]
                       for t in CHLC_GRID], dtype=float)
    roots = []
    for i in range(len(CHLC_GRID) - 1):
        a, b = float(vals[i]), float(vals[i + 1])
        lo, hi = float(CHLC_GRID[i]), float(CHLC_GRID[i + 1])
        if abs(a) < 1e-10:
            roots.append(lo)
        elif a * b < 0:
            roots.append(float(brentq(
                lambda t: evaluate(ctx, t, f_c, shortwave)[0],
                lo, hi, xtol=1e-12
            )))
    if abs(vals[-1]) < 1e-10:
        roots.append(float(CHLC_GRID[-1]))
    roots = sorted(roots)
    unique = []
    for root in roots:
        if not unique or abs(root - unique[-1]) > 1e-8:
            unique.append(root)
    # Also collect PAR for NO_ROOT / MULTIPLE_ROOTS, using a diagnostic theta.
    # This does not invent a balanced Chl:C where no unique solution exists.
    if len(unique) != 1:
        _, diagnostic_par = evaluate(ctx, 0.025, f_c, shortwave)
        return (float('nan'), diagnostic_par, float('nan'),
                'NO_ROOT' if not unique else 'MULTIPLE_ROOTS')
    theta = float(unique[0])
    residual, par_wm2 = evaluate(ctx, theta, f_c, shortwave)
    if abs(residual) > BALANCE_TOL_PER_DAY:
        raise RuntimeError(f'Non-zero balance residual: {residual} d^-1'
                           f' at fC={f_c}, SW={shortwave}')
    return theta, par_wm2, residual, 'ONE_ROOT'


def validate_reference_against_previous(rows):
    """Optional audit: ensure f_C=0.24 agrees with passed Legacy/ref CSV."""
    if not VERIFIED_REFERENCE.is_file():
        print('Reference cross-check: prior consistency CSV absent (skipping).')
        return
    with VERIFIED_REFERENCE.open(newline='') as f:
        prev = list(csv.DictReader(f))
    reference = [r for r in rows if r['fC'] == FC_REF]
    if len(prev) != len(reference):
        raise RuntimeError('Previous reference CSV length differs from this experiment')
    deltas = []
    for current, earlier in zip(reference, prev):
        if not np.isclose(current['shortwave_W_m2'],
                          float(earlier['shortwave_W_m2']), rtol=0, atol=1e-9):
            raise RuntimeError('Reference cross-check shortwave grids differ')
        if current['status'] != 'ONE_ROOT':
            raise RuntimeError('Reference Chl:C root missing; prior audit found one')
        deltas.append(abs(current['balanced_ChlC_g_g'] -
                          float(earlier['reference_allocation_balanced_ChlC_g_g'])))
    max_delta = max(deltas)
    print(f'Previous verified-reference cross-check: max |Delta Chl:C|={max_delta:.3e}')
    if max_delta > 1e-7:
        raise RuntimeError('New reference curve differs from the previously verified curve')


def main():
    if not YAML.is_file():
        raise FileNotFoundError(YAML)
    if FC_REF not in FC_VALUES:
        raise ValueError('Reference f_C must be included')
    if not all(0.05 <= fc <= 0.55 and 1.0-FLH_FIXED-fc > 0
               for fc in FC_VALUES):
        raise ValueError('Prescribed f_C violates allocation bounds/budget')

    x_obs, y_obs = read_observations()
    ctx = build_model()
    rows = []
    for fc in FC_VALUES:
        for shortwave in SW_GRID:
            theta, par_w, residual, status = solve_balanced_chlc(
                ctx, float(fc), float(shortwave))
            rows.append({
                'fC': float(fc),
                'fV_implied': float(1.0-FLH_FIXED-fc),
                'shortwave_W_m2': float(shortwave),
                'model_PAR_W_m2': par_w,
                'model_PAR_umol_photons_m2_s': par_w * PAR_UMOL_PER_J,
                'balanced_ChlC_g_g': theta,
                'balance_residual_d_inv': residual,
                'status': status,
            })
        cases = [r for r in rows if r['fC'] == fc]
        n = sum(r['status'] == 'ONE_ROOT' for r in cases)
        print(f'fC={fc:.2f} (implied fV={1-FLH_FIXED-fc:.2f}): '
              f'{n}/{len(SW_GRID)} unique balanced-growth roots', flush=True)
        if n < len(SW_GRID):
            print(f'  WARNING: {len(SW_GRID)-n} conditions lack one unique root. '
                  'These will appear as gaps, NOT connected or interpolated.')
    validate_reference_against_previous(rows)

    reference = {r['shortwave_W_m2']: r for r in rows if r['fC'] == FC_REF}
    for fc in FC_VALUES:
        if fc == FC_REF:
            continue
        delta = [abs(r['balanced_ChlC_g_g'] -
                     reference[r['shortwave_W_m2']]['balanced_ChlC_g_g'])
                 for r in rows if r['fC'] == fc and r['status'] == 'ONE_ROOT'
                 and reference[r['shortwave_W_m2']]['status'] == 'ONE_ROOT']
        if delta:
            print(f'Against fC={FC_REF:.2f}, fC={fc:.2f}: '
                  f'max |Delta Chl:C| at matched shortwave={max(delta):.5g} g/g; '
                  f'{len(delta)} paired conditions')

    OUTPUT.mkdir(parents=True, exist_ok=True)
    csv_out = OUTPUT / 'ersem_alloc_fixed_fc_sensitivity.csv'
    with csv_out.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, ax = plt.subplots(figsize=(8.2, 5.6), constrained_layout=True)
    colors = {0.12: '#327B97', 0.24: '#2255A4',
              0.38: '#CF7833', 0.50: '#8B4F90'}
    for fc in FC_VALUES:
        case = [r for r in rows if r['fC'] == fc]
        # Keep the original shortwave order; NaN gaps prevent a deceptive
        # curve connecting separately valid branches through invalid points.
        x = np.asarray([r['model_PAR_umol_photons_m2_s'] for r in case])
        y = np.asarray([r['balanced_ChlC_g_g'] for r in case])
        ax.plot(x, y, lw=(2.8 if fc == FC_REF else 1.9),
                color=colors[fc], ls=('-' if fc == FC_REF else '--'),
                label=rf'Fixed $f_C={fc:.2f}$' +
                      (' (Legacy-equivalent reference)' if fc == FC_REF else ''))
    ax.scatter(x_obs, y_obs, color='black', s=42, zorder=5, marker='o',
               label=r'$T.$ pseudonana (Geider et al., 1997; digitised)')
    finite_chlc = np.asarray([r['balanced_ChlC_g_g'] for r in rows])
    upper = max(0.065, 1.07 * float(np.nanmax(np.r_[finite_chlc, y_obs])))
    ax.set(xlim=(0, 600), ylim=(0, upper),
           xlabel=r'Irradiance ($\mu$mol photons m$^{-2}$ s$^{-1}$)',
           ylabel=r'Chlorophyll $a$:carbon (g Chl $a$ g$^{-1}$ C)')
    ax.grid(alpha=0.15)
    ax.legend(frameon=False, fontsize=9, loc='upper right')
    ax.text(0.018, 0.024,
            'P1 diatoms; 291 K; fixed initial quotas; prescribed $f_C$ (not optimised)',
            transform=ax.transAxes, fontsize=8.5, color='0.32')
    png = OUTPUT / 'ersem_alloc_fixed_fc_sensitivity.png'
    pdf = OUTPUT / 'ersem_alloc_fixed_fc_sensitivity.pdf'
    fig.savefig(png, dpi=250)
    fig.savefig(pdf)
    plt.close(fig)

    print('\n===== CONDITIONAL fC SENSITIVITY COMPLETE =====')
    print('P1 generic diatoms; 291 K; same initial states at each evaluation.')
    print('Only fC and chlorophyll varied; all other model states fixed.')
    print('Approximate PAR conversion: 4.57 umol photons/J.')
    print('NOT a nutrient-quota equilibrium or optimised/prognostic fC solution.')
    print('CSV   :', csv_out)
    print('Figure:', png)
    print('PDF   :', pdf)


if __name__ == '__main__':
    main()
