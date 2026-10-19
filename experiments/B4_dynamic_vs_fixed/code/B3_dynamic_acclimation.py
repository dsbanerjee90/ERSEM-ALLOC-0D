# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
#
# Experiment B3 — time-dependent N/P/PAR acclimation of ERSEM-ALLOC P1.
#
# IMPORTANT:
# - A conservative, explicit-Euler, standalone PyFABM test harness, NOT a
#   reproduction of the NEMO host time integration or a steady-state solver.
# - Evolve all FABM prognostic states, plus ONLY P1's host-owned fC dependency.
#   P2-P4 allocations remain fixed at the Legacy-equivalent reference 0.24.
# - External NO3, NH4, PO4 and Si remain prescribed CONSTANT BATH values via
#   state resetting after each substep (open chemostat-type nutrient reservoirs).
#   Nutrient budgets therefore include implicit external sources/sinks.
# - Temperature, irradiance and other environmental dependencies fixed.
# - Never edit source, YAML, original experiment dirs, or earlier outputs.
# - Stop on negative/new nonfinite concentrations: NO silent clipping.
# - EXCEPTION: G2/o_deep is a signed deep benthic oxygen profile integral,
#   explicitly permitted by ERSEM benthic_column_dissolved_matter.F90.
#   No other state is exempted by this experiment.
# - 'Quasi-stationary' heuristic is NOT proof of stability or optimality.
# - Run --days 0.5 first, then --days 7 --dt 300, and repeat at --dt 150.

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from section2_pyfabm_init import build_model

SEC_DAY = 86400.0
F0 = 0.24
FC_MIN = 0.05
FC_MAX = 0.55
FLH = 0.38
# Fixed nutrient bath values are native ERSEM mmol N/P/Si m^-3 concentrations.
# These are experimental contrasts, not proven physiological classifications.
TREATMENTS = {
    'replete': (8.0, 0.4),
    'low_N': (0.03, 0.4),
    'low_P': (8.0, 0.003),
    'low_NP': (0.03, 0.003),
}
NH4_BATH = 0.001  # prevent an unchanged NH4 pool from rescuing N-poor treatments
SI_BATH = 4.5     # hold silicate equal in all treatments
QNL_MIN = 0.0042  # confirmed P1/qnlc, mmol N / mg C
QPL_MIN = 0.0001  # confirmed P1/qplc, mmol P / mg C
PAR_UMOL_PER_J = 4.57  # approximate, spectrum-dependent conversion
REQUIRED_DIAG = ('G_C', 'G_Q', 'Gamma_C', 'dfCdt', 'condJ')
OPTIONAL_DIAG = ('fO3PIc', 'fN3PIn', 'fN4PIn', 'fN1PIp')
STATE_NAMES = ('P1/c', 'P1/n', 'P1/p', 'P1/Chl', 'N3/n', 'N4/n', 'N1/p', 'N5/s')
# Verified source: ERSEM's benthic_column_dissolved_matter.F90 uses
# `if (.not.info%nonnegative)` when registering/integrating the deep integral.
# The initial G2/o_deep=0 and its negative tendency are intentional in this branch.
# Do NOT generalise this exception to all *_deep states.
ALLOWED_SIGNED_STATES = frozenset(('G2/o_deep',))
C_COLORS = {'replete': '#252f82', 'low_N': '#cf7936',
            'low_P': '#187d89', 'low_NP': '#913d83'}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--days', type=float, default=0.5,
                   help='duration in days; default 0.5 for safe pilot')
    p.add_argument('--dt', type=float, default=300.,
                   help='explicit integration timestep in seconds')
    p.add_argument('--lights', default='20,100',
                   help='shortwave fluxes in W/m2, comma separated')
    p.add_argument('--output', type=Path, default=HERE.parent / 'results' / 'np_par_acclimation_B3',
                   help='NEW output directory, must not already contain results')
    p.add_argument('--save-hours', type=float, default=1.0)
    args = p.parse_args()
    if args.days <= 0 or args.dt <= 0 or args.save_hours <= 0:
        p.error('days, dt and save-hours must be positive')
    nstep = args.days * SEC_DAY / args.dt
    if not math.isclose(nstep, round(nstep), abs_tol=1e-8):
        p.error('days*86400/dt must be an integer: choose matching duration/timestep')
    nsave = args.save_hours * 3600 / args.dt
    if not math.isclose(nsave, round(nsave), abs_tol=1e-8):
        p.error('save-hours*3600/dt must be an integer')
    try:
        args.lights = tuple(float(s) for s in args.lights.split(','))
    except ValueError:
        p.error('lights must be comma-separated numbers, e.g. 20,100')
    if not args.lights or any(x < 0 for x in args.lights):
        p.error('invalid shortwave list')
    args.nsteps = int(round(nstep))
    args.save_every = int(round(nsave))
    return args


def diagnose(m, names):
    answer = {}
    for name in names:
        v = m.findDiagnosticVariable('P1/' + name)
        answer[name] = v
    return answer


def dump_csv(path, rows):
    if not rows:
        return
    if path.exists():
        raise FileExistsError(f'Will not replace previous result: {path}')
    with path.open('w', newline='') as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def sample(m, state, fc, sw, k, dt, diags, parvar, treatment, bath,
           nidx, fdep, par_conversion):
    m.state[:] = state
    fdep.value = fc
    rates = np.asarray(m.getRates(), dtype=float).copy()
    if rates.shape != state.shape or not np.isfinite(rates).all():
        raise RuntimeError(f'FABM rates invalid at {treatment}, step {k}')
    dg = {key: float(var.value) if var is not None else float('nan')
          for key, var in diags.items()}
    if not all(np.isfinite(dg[key]) for key in REQUIRED_DIAG):
        raise RuntimeError(f'Missing/invalid allocation gradients at {treatment}, step {k}')
    if not np.isclose(dg['Gamma_C'], dg['G_C']+dg['G_Q'],
                      rtol=1e-7, atol=1e-6):
        raise RuntimeError(f'Gamma != G_C+G_Q at {treatment}, step {k}')
    if abs(dg['dfCdt']) > 1.500001:
        raise RuntimeError('Allocation rate exceeds configured 1.5 day^-1 limit')
    c = state[nidx['P1/c']]
    n = state[nidx['P1/n']]
    p = state[nidx['P1/p']]
    if min(c, n, p) <= 0:
        raise RuntimeError(f'Zero or negative P1 C/N/P at {treatment}, step {k}')
    phot = float(parvar.value)
    if not np.isfinite(phot) or phot < 0:
        raise RuntimeError('Invalid PAR diagnostic')
    row = {
        't_days': k * dt / SEC_DAY,
        'treatment': treatment, 'shortwave_W_m2': sw,
        'PAR_W_m2': phot, 'PAR_umol_m2_s_approx': phot*par_conversion,
        'bath_NO3': bath['N3/n'], 'bath_NH4': bath['N4/n'],
        'bath_PO4': bath['N1/p'], 'bath_Si': bath['N5/s'],
        'fC': fc, 'fV': 1.-FLH-fc,
        'P1_C': c, 'P1_N': n, 'P1_P': p,
        'P1_Chl': state[nidx['P1/Chl']],
        'N_to_C_native': n/c, 'P_to_C_native': p/c,
        'below_native_N_quota_min': int(n/c < QNL_MIN),
        'below_native_P_quota_min': int(p/c < QPL_MIN),
        'Chl_to_C_native': state[nidx['P1/Chl']]/c,
        'G_C': dg['G_C'], 'G_Q': dg['G_Q'],
        'Gamma_C': dg['Gamma_C'], 'dfCdt_d_inv': dg['dfCdt'],
        'condJ': dg['condJ'],
        'P1_gross_C_fixation': dg['fO3PIc'],
        'P1_NO3_uptake': dg['fN3PIn'],
        'P1_NH4_uptake': dg['fN4PIn'],
        'P1_PO4_uptake': dg['fN1PIp'],
        'gross_C_fixation_per_C': dg['fO3PIc']/c,
        'near_fC_bound': int(abs(fc-FC_MIN)<1e-8 or abs(fc-FC_MAX)<1e-8),
    }
    return row, rates


def run_case(args, sw, treatment, bathN, bathP):
    ctx = build_model()
    m, original, _, _, fdep, shortwave, parvar = ctx
    nidx = {v.name: i for i, v in enumerate(m.state_variables)}
    missing = [name for name in STATE_NAMES if name not in nidx]
    if missing:
        raise RuntimeError(f'Expected model states absent: {missing}')
    diags = diagnose(m, REQUIRED_DIAG + OPTIONAL_DIAG)
    missing_diag = [name for name in REQUIRED_DIAG if diags[name] is None]
    if missing_diag:
        raise RuntimeError(f'Missing required diagnostics: {missing_diag}')
    # Clamp the specified external dissolved pools. This is OPEN to an
    # implicit reservoir, not a mass-conserving closed-system experiment.
    bath = {'N3/n': bathN, 'N4/n': NH4_BATH,
            'N1/p': bathP, 'N5/s': SI_BATH}
    for name, value in bath.items():
        if value < 0:
            raise ValueError(f'Negative bath {name}')
    bath_indices = np.array([nidx[name] for name in bath], dtype=int)
    bath_values = np.array(list(bath.values()), dtype=float)
    state = original.copy()
    state[bath_indices] = bath_values
    if not np.isfinite(state).all():
        raise RuntimeError('Nonfinite starting state')
    # Default conservative protection: forbid negative values for states that
    # begin >= 0, except source-verified signed states. The exception is NOT
    # a numerical clip; its Euler-computed negative value remains unchanged.
    absent_signed = ALLOWED_SIGNED_STATES.difference(nidx)
    if absent_signed:
        raise RuntimeError(f'Expected signed benthic states absent: {sorted(absent_signed)}')
    signed_indices = np.array([nidx[name] for name in sorted(ALLOWED_SIGNED_STATES)], dtype=int)
    nonneg_indices = np.array([i for i in np.flatnonzero(state >= 0.)
                               if i not in set(signed_indices)], dtype=int)
    min_signed = state[signed_indices].copy()
    fc = F0
    shortwave.value = float(sw)
    rows = []
    fc_bounds_hits = 0
    min_cond = float('inf')
    # Host-owned fC has a 1/day tendency; FABM state tendencies are per second.
    for k in range(args.nsteps+1):
        row, rates = sample(m, state, fc, sw, k, args.dt, diags, parvar,
                            treatment, bath, nidx, fdep, PAR_UMOL_PER_J)
        min_cond = min(min_cond, row['condJ'])
        if k % args.save_every == 0 or k == args.nsteps:
            rows.append(row)
        if k == args.nsteps:
            break
        nxt = state + args.dt*rates
        nxt[bath_indices] = bath_values
        min_signed = np.minimum(min_signed, nxt[signed_indices])
        if not np.isfinite(nxt).all():
            raise RuntimeError(f'Nonfinite state after t={(k+1)*args.dt/SEC_DAY:.5g}d '
                               f'in {treatment}, SW={sw}; try --dt 150 or 30')
        negatives = nonneg_indices[nxt[nonneg_indices] < 0.]
        if negatives.size:
            j = negatives[np.argmin(nxt[negatives])]
            state_name = m.state_variables[int(j)].name
            raise RuntimeError(f'Negative state {state_name}={nxt[j]:.5g} at '
                               f't={(k+1)*args.dt/SEC_DAY:.5g}d in '
                               f'{treatment}, SW={sw}. No clipping applied. '
                               'Inspect the state definition and its biological '
                               'rate before changing dt or adding exceptions.')
        # Reproduce the model's bounded, rate-limited host-owned P1 allocation
        # as a forward Euler dependency update. Diagnostic already rate-capped
        # to +/-1.5 day^-1 in native Fortran.
        raw_next_fc = fc + (args.dt/SEC_DAY)*row['dfCdt_d_inv']
        new_fc = float(np.clip(raw_next_fc, FC_MIN, FC_MAX))
        fc_bounds_hits += int(new_fc != raw_next_fc)
        state, fc = nxt, new_fc
    return rows, fc_bounds_hits, min_cond, dict(zip(sorted(ALLOWED_SIGNED_STATES),
                                                     map(float, min_signed)))


def plot_light(out, sw, by_treatment):
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.2), constrained_layout=True)
    keys = [('fC', r'Carbon machinery fraction $f_C$'),
            ('N_to_C_native', r'P1 intracellular N:C (native units)'),
            ('P_to_C_native', r'P1 intracellular P:C (native units)'),
            ('Gamma_C', r'Allocation gradient $G_C+G_Q$')]
    for ax, (key, label) in zip(axes.flat, keys):
        for tr, rows in by_treatment.items():
            x = [r['t_days'] for r in rows]
            y = [r[key] for r in rows]
            ax.plot(x, y, color=C_COLORS[tr], lw=1.8, label=tr.replace('_', ' '))
        ax.set(ylabel=label, xlabel='Elapsed time (days)')
        ax.grid(alpha=.18)
        if key == 'fC':
            ax.axhline(F0, ls=':', color='0.4', lw=1)
            ax.axhline(FC_MIN, ls=':', color='0.8', lw=.8)
            ax.axhline(FC_MAX, ls=':', color='0.8', lw=.8)
        if key == 'Gamma_C':
            ax.axhline(0, ls=':', color='0.4', lw=1)
    axes[0, 0].legend(loc='best', frameon=False, fontsize=9)
    fig.suptitle(f'ERSEM-ALLOC B3 — dynamic P1 allocation under fixed nutrient baths; '
                 f'SW={sw:g} W m$^{{-2}}$', fontsize=13)
    fig.text(.5, -.015,
             'Full PyFABM state evolves; NO3, NH4, PO4 and Si externally prescribed. '
             'Constant light. Trajectories are not proven physiological optima.',
             fontsize=8, ha='center')
    stem = f'np_acclimation_SW{sw:g}'.replace('.', 'p')
    png, pdf = out / f'{stem}.png', out / f'{stem}.pdf'
    for p in (png, pdf):
        if p.exists():
            raise FileExistsError(f'Will not replace {p}')
    fig.savefig(png, dpi=260, bbox_inches='tight')
    fig.savefig(pdf, bbox_inches='tight')
    plt.close(fig)
    return png, pdf


def main():
    args = parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise FileExistsError(f'Output directory contains existing files: {args.output}; '
                              'choose --output <new-directory>')
    args.output.mkdir(parents=True, exist_ok=True)
    print('===== ERSEM-ALLOC B3 TIME-DEPENDENT N/P/PAR ACCLIMATION =====', flush=True)
    print(f'Duration {args.days} days; Euler dt={args.dt:g}s; '
          f'{args.nsteps} steps / treatment; SW={args.lights}', flush=True)
    print('Full FABM state integrated except external nitrate, ammonium, phosphate, silicate.')
    print('Only P1 fC evolves; P2..P4 fC remain at reference=0.24.')
    print('N and P treatments are ambient BATH contrasts, not verified quota-limitation labels.')
    print('Source / YAML untouched; original 3D outputs untouched.')
    print('Source-verified signed-state exception: G2/o_deep (deep benthic oxygen integral).')
    print('All other originally nonnegative states retain the stop-on-negative check.')
    summaries = []
    for sw in args.lights:
        by_tr = {}
        for tr, (n, p) in TREATMENTS.items():
            print(f'RUN sw={sw:g}, tr={tr}, NO3={n:g}, PO4={p:g}', flush=True)
            rows, fc_hits, min_cond, signed_min = run_case(args, sw, tr, n, p)
            sw_label = f'{sw:g}'.replace('.', 'p')
            csv_name = f'np_acclimation_SW{sw_label}_{tr}.csv'
            csv_path = args.output / csv_name
            dump_csv(csv_path, rows)
            by_tr[tr] = rows
            tail = [r for r in rows if r['t_days'] >= 0.75*args.days]
            change = abs(tail[-1]['fC']-tail[0]['fC']) if len(tail)>1 else float('nan')
            summaries.append({
                'shortwave_W_m2': sw, 'treatment': tr,
                'bath_NO3': n, 'bath_NH4': NH4_BATH, 'bath_PO4': p,
                'days': args.days, 'dt_seconds': args.dt,
                'fC_initial': F0, 'fC_final': rows[-1]['fC'],
                'fC_change_last_25pct_days': change,
                'N_to_C_final': rows[-1]['N_to_C_native'],
                'P_to_C_final': rows[-1]['P_to_C_native'],
                'Gamma_final': rows[-1]['Gamma_C'],
                'G_C_final': rows[-1]['G_C'], 'G_Q_final': rows[-1]['G_Q'],
                'min_condJ': min_cond, 'fC_bounds_clamps': fc_hits,
                'min_G2_o_deep_signed': signed_min['G2/o_deep'],
                'P1_C_final': rows[-1]['P1_C'],
                'saved_samples_below_N_quota_min': sum(r['below_native_N_quota_min'] for r in rows),
                'saved_samples_below_P_quota_min': sum(r['below_native_P_quota_min'] for r in rows),
                'P1_fixation_final': rows[-1]['P1_gross_C_fixation'],
                'status': 'INTEGRATED_NOT_VALIDATED_STABLE',
            })
            print(f"  end fC={rows[-1]['fC']:.5f}, qN={rows[-1]['N_to_C_native']:.5g}, "
                  f"qP={rows[-1]['P_to_C_native']:.5g}, "
                  f"Gamma={rows[-1]['Gamma_C']:+.4g}, min condJ={min_cond:.3g}, "
                  f'boundary clamp events={fc_hits}; '
                  f"signed G2/o_deep min={signed_min['G2/o_deep']:+.4g}; "
                  f'low-quota samples N/P='
                  f"{sum(r['below_native_N_quota_min'] for r in rows)}/"
                  f"{sum(r['below_native_P_quota_min'] for r in rows)}", flush=True)
        png, pdf = plot_light(args.output, sw, by_tr)
        print('FIGURES:', png, pdf, flush=True)
    dump_csv(args.output / 'np_acclimation_B3_summary.csv', summaries)
    print('SUMMARY:', args.output / 'np_acclimation_B3_summary.csv')
    print('Do not call endpoint optimum or stable without convergence, branch and dt tests.')


if __name__ == '__main__':
    main()
