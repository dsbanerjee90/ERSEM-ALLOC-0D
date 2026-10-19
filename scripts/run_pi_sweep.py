# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
#
# Controlled PyFABM photosynthesis-irradiance (P-I) sensitivity experiment:
# P1 gross inorganic-carbon fixation vs PAR for prescribed machinery f_C.
#
# Two genuinely distinct tests:
#  A. Same P1 carbon, same Chl:C, same other states (isolates the fixed
#     carbon-allocation effect without a chlorophyll-acclimation confound).
#  B. Same P1 carbon and other states, but set Chl:C to its previously
#     calculated balanced-growth value separately at each f_C and PAR.
#
# The GPP proxy is NOT net d(P1/c)/dt. The exact gross inorganic-C -> P1
# carbon flux diagnostic P1/fO3PIc is used. If it is unavailable, STOP
# rather than silently substitute net growth or a guessed process.
#
# This experiment does not infer nutrient starvation, solve optimal f_C,
# or estimate ERSEM's intrinsic alpha parameter. No existing files or
# model configurations are changed; outputs are all new files.

from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyfabm

ROOT = Path(__file__).resolve().parents[1]  # This new PAPER workspace
YAML = ROOT / 'configs/fabm.yaml'
OUT = ROOT / 'results'
BALANCED_CSV = OUT / 'ersem_alloc_fixed_fc_sensitivity.csv'

FC_REF = 0.24
FC_VALUES = (0.12, 0.24, 0.38, 0.50)
FIXED_CHLC = 0.025  # g Chl a per g C for ALL f_C and light levels, panel A.
PAR_UMOL_PER_J = 4.57  # approximate spectral conversion; as in earlier test.
TEMP_K = 291.0

# Include low-irradiance points to examine initial P-I behaviour; do not
# claim these slopes are the ERSEM model parameter alpha.
SW_FIXED = np.unique(np.r_[np.linspace(0., 14., 29),
                             np.linspace(16., 350., 64)])

ENV = {
    'cell_thickness': 1.0,
    'temperature': TEMP_K - 273.15,
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


def make_model():
    if not YAML.is_file():
        raise FileNotFoundError(f'Original read-only FABM YAML not found: {YAML}')
    m = pyfabm.Model(str(YAML))
    missing = sorted({d.name for d in m.dependencies} - set(ENV))
    if missing:
        raise RuntimeError(f'Undefined FABM dependencies: {missing}')
    for d in m.dependencies:
        d.value = ENV[d.name]
    m.setCellThickness(1.0)
    if not m.checkReady():
        raise RuntimeError('FABM dependency check failed')

    state_names = [v.name for v in m.state_variables]
    carbon_idx, chl_idx = state_names.index('P1/c'), state_names.index('P1/Chl')
    initial = np.asarray(m.state, dtype=float).copy()
    if not np.isfinite(initial).all() or initial[carbon_idx] <= 0:
        raise RuntimeError('Nonfinite or nonpositive initial P1 carbon')

    fc = m.findDependency('P1/fC_state')
    sw = m.findDependency('surface_downwelling_shortwave_flux')
    par = m.findDiagnosticVariable('light/parEIR')
    gpp = m.findDiagnosticVariable('P1/fO3PIc')
    if gpp is None:
        matching = [v.name for v in m.diagnostic_variables
                    if v.name.startswith('P1/') and any(
                        s in v.name.lower() for s in
                        ('fo3', 'pp', 'photo', 'fix', 'gpp', 'cfix', 'assim', 'gross'))]
        raise RuntimeError(
            'Cannot find ERSEM P1 gross inorganic-carbon fixation diagnostic '
            "'P1/fO3PIc'. No net-tracer-rate substitution was made.\n"
            f'Candidate P1 diagnostic names: {matching}\n'
            'Please send the printed candidates and the P1 flux registration '
            'from primary_producer.F90 for a precise adjustment.'
        )
    if any(v is None for v in (fc, sw, par)):
        raise RuntimeError('Missing P1/fC_state, shortwave, or light/parEIR')
    units = str(getattr(gpp, 'units', '') or 'model-native flux units')
    long_name = str(getattr(gpp, 'long_name', '') or '')
    print(f'Using ERSEM diagnostic: {gpp.name}; units={units}; description={long_name}',
          flush=True)
    print(f'Initial P1/c = {initial[carbon_idx]:.9g}, Chl:C held at '
          f'{FIXED_CHLC} in panel A.', flush=True)
    return m, initial, carbon_idx, chl_idx, fc, sw, par, gpp, units


def evaluate(ctx, fc_value, shortwave, chlc):
    m, initial, ic, ih, fc, sw, par, gpp, _ = ctx
    # Reset full state on every call: no integration, no state history.
    m.state[:] = initial
    m.state[ih] = float(initial[ic]) * float(chlc)
    fc.value = float(fc_value)
    sw.value = float(shortwave)
    rates = np.asarray(m.getRates(), dtype=float)
    if not np.isfinite(rates).all():
        raise RuntimeError('Nonfinite FABM state rates')
    gpp_value = float(gpp.value)
    par_w = float(par.value)
    if not np.isfinite(gpp_value) or not np.isfinite(par_w) or par_w < 0:
        raise RuntimeError(f'Nonfinite ERSEM GPP or PAR: fc={fc_value}, sw={shortwave}')
    return {
        'fC': float(fc_value),
        'fV_implied': float(1.0 - 0.38 - fc_value),
        'shortwave_W_m2': float(shortwave),
        'PAR_W_m2': par_w,
        'PAR_umol_photons_m2_s': par_w * PAR_UMOL_PER_J,
        'ChlC_g_g': float(chlc),
        'P1_carbon_native_concentration': float(m.state[ic]),
        'gross_C_fixation_P1_fO3PIc_native': gpp_value,
        'gross_C_fixation_per_P1C_native': gpp_value / float(m.state[ic]),
    }


def load_balanced_cases():
    if not BALANCED_CSV.is_file():
        print('Previous balanced Chl:C CSV not found; panel B will be omitted.',
              flush=True)
        return []
    with BALANCED_CSV.open(newline='') as f:
        rows = list(csv.DictReader(f))
    cases = []
    for r in rows:
        if r.get('status') != 'ONE_ROOT':
            continue
        fc = float(r['fC'])
        chl = float(r['balanced_ChlC_g_g'])
        sw = float(r['shortwave_W_m2'])
        if fc not in FC_VALUES or not (np.isfinite(chl) and chl > 0):
            continue
        cases.append((fc, sw, chl))
    print(f'Balanced-growth cases read: {len(cases)} from {BALANCED_CSV.name}',
          flush=True)
    return cases


def add_curves(ax, rows, field, label_max=False):
    colors = {0.12: '#327B97', 0.24: '#2255A4',
              0.38: '#CF7833', 0.50: '#8B4F90'}
    for fc in FC_VALUES:
        case = sorted((r for r in rows if r['fC'] == fc),
                      key=lambda r: r['PAR_umol_photons_m2_s'])
        if not case:
            continue
        label = rf'$f_C={fc:.2f}$'
        if fc == FC_REF:
            label += ' (reference)'
        ax.plot([r['PAR_umol_photons_m2_s'] for r in case],
                [r[field] for r in case],
                color=colors[fc], lw=(2.7 if fc == FC_REF else 2.1),
                ls=('-' if fc == FC_REF else '--'), label=label)
    ax.set_xlabel(r'PAR ($\mu$mol photons m$^{-2}$ s$^{-1}$)')
    ax.grid(alpha=0.18)
    ax.legend(frameon=False, fontsize=9)
    ax.set_xlim(left=0, right=600)


def print_diagnostics(fixed, balanced):
    print('\n===== P1 CARBON-FIXATION P-I SENSITIVITY =====', flush=True)
    print('A: identical carbon, chlorophyll, background N/P, temperature, light;',
          'only fC prescribed.', flush=True)
    for fc in FC_VALUES:
        case = sorted([r for r in fixed if r['fC'] == fc],
                      key=lambda r: r['PAR_umol_photons_m2_s'])
        if not case:
            continue
        x = np.asarray([r['PAR_umol_photons_m2_s'] for r in case])
        y = np.asarray([r['gross_C_fixation_P1_fO3PIc_native'] for r in case])
        # Low-light secant estimate, not a claim about the parameter alpha.
        inds = (x >= 0) & (x <= 20)
        if inds.sum() >= 4:
            xx = x[inds] - x[inds][0]
            yy = y[inds] - y[inds][0]
            low_light_slope = np.dot(xx, yy) / np.dot(xx, xx)
        else:
            low_light_slope = float('nan')
        y300 = np.interp(300., x, y)
        print(f'  fC={fc:.2f}: GPP(PAR=300)={y300:.8e} model-native; '
              f'low-PAR secant slope={low_light_slope:.8e} '
              '(native flux per umol photons)', flush=True)
    if balanced:
        print('B: Chl:C separately balanced at each irradiance (from previous '
              'sensitivity CSV); all other states fixed.', flush=True)
    print('NOTE: P1/fO3PIc is gross inorganic-C-to-diatom-C production, NOT '
          'net P1 carbon growth.', flush=True)
    print('NOTE: fC is prescribed, NOT an optimiser solution or N-starvation '
          'experiment.', flush=True)
    print('NOTE: 4.57 umol photons/J PAR conversion is approximate.', flush=True)
    print('NOTE: Low-PAR fitted slopes are diagnostic response slopes, NOT '
          'necessarily the ERSEM photosynthetic efficiency parameter alpha.', flush=True)


def main():
    ctx = make_model()
    native_units = ctx[-1]
    if not all(0.05 <= f <= 0.55 for f in FC_VALUES):
        raise ValueError('fC outside ERSEM-ALLOC prescribed bounds')
    fixed = []
    for fc in FC_VALUES:
        for sw in SW_FIXED:
            fixed.append(evaluate(ctx, fc, sw, FIXED_CHLC))
        print(f'fC={fc:.2f}: {len(SW_FIXED)} fixed-Chl:C PAR points', flush=True)
    balanced = []
    for fc, sw, chl in load_balanced_cases():
        balanced.append(evaluate(ctx, fc, sw, chl))
    print_diagnostics(fixed, balanced)

    OUT.mkdir(parents=True, exist_ok=True)
    csvpath = OUT / 'ersem_alloc_p1_photosynthesis_irradiance.csv'
    combined = [dict(mode='FIXED_CHLC', **r) for r in fixed] + [
        dict(mode='BALANCED_CHLC', **r) for r in balanced]
    with csvpath.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(combined[0].keys()))
        writer.writeheader()
        writer.writerows(combined)

    fig, axes = plt.subplots(1, 2 if balanced else 1, figsize=(12.4, 5.2)
                              if balanced else (7.4, 5.2),
                              constrained_layout=True)
    if not balanced:
        axes = [axes]
    add_curves(axes[0], fixed, 'gross_C_fixation_P1_fO3PIc_native')
    axes[0].set_title('A  Fixed Chl:C = 0.025; identical biomass and quotas'
                      if balanced else 'P1 carbon fixation vs irradiance')
    axes[0].set_ylabel(f'P1 gross carbon-fixation flux ({native_units})')
    if balanced:
        add_curves(axes[1], balanced, 'gross_C_fixation_P1_fO3PIc_native')
        axes[1].set_title('B  Chl:C balanced separately at each irradiance')
        axes[1].set_ylabel(f'P1 gross carbon-fixation flux ({native_units})')
    fig.suptitle('ERSEM-ALLOC P1: prescribed machinery allocation and P–I response',
                 fontsize=13)
    fig.text(.5, .012,
             '291 K; constant initial N/P states; same P1 carbon; fC prescribed,',
             ha='center', fontsize=8.5, color='0.35')
    png = OUT / 'ersem_alloc_p1_photosynthesis_irradiance.png'
    pdf = OUT / 'ersem_alloc_p1_photosynthesis_irradiance.pdf'
    fig.savefig(png, dpi=250)
    fig.savefig(pdf)
    plt.close(fig)
    print('\nCSV   :', csvpath, '\nPNG   :', png, '\nPDF   :', pdf, flush=True)


if __name__ == '__main__':
    main()
