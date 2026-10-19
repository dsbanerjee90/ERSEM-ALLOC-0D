# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
#
# ERSEM-ALLOC Experiment B4: matched 0D dynamic-fC vs fixed-reference-fC controls.
# Requires B3_dynamic_acclimation.py and
# section2_pyfabm_init.py in the SAME scripts directory.
#
# This is a diagnostic, fixed-nutrient-bath, full-FABM Euler integration test.
# The same coupled biogeochemical states evolve in each arm. Only the host-owned
# P1/fC_state is updated in the dynamic arm; the fixed arm holds it at 0.24.
# The fC capacity multipliers remain enabled in BOTH arms. This is NOT a
# separately compiled Legacy run, a conservation-closed system, a proven
# physiological optimum, or an exact reproduction of the NEMO host integrator.
# Benthic state G2/o_deep is source-confirmed to be signed; it is NOT clipped.
# All existing source, YAML, and experimental output files are left untouched.

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import B3_dynamic_acclimation as b3

QN_REDFIELD = 0.0126  # ersem/shared.F90: qnRPIcX, mmol N/mg C
QP_REDFIELD = 0.000786  # ersem/shared.F90: qpRPIcX, mmol P/mg C
SELECTED_TREATMENTS = ('low_N', 'low_P', 'low_NP')
SOURCE_VERIFIED_SIGNED = frozenset(('G2/o_deep',))
NEEDED_FLUXES = ('fO3PIc', 'fN3PIn', 'fN4PIn', 'fN1PIp')


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--days', type=float, default=0.5)
    ap.add_argument('--dt', type=float, default=150.0,
                    help='Euler timestep seconds; use 150 for the main run')
    ap.add_argument('--lights', default='20,100',
                    help='comma-separated shortwave W/m2')
    ap.add_argument('--save-hours', type=float, default=1.0)
    ap.add_argument('--include-replete', action='store_true',
                    help='Optional separate robustness check; known gradient switching')
    ap.add_argument('--output', type=Path,
                    default=b3.HERE.parent / 'results' / 'B4_dynamic_vs_fixed_pilot')
    a = ap.parse_args()
    if not all(math.isfinite(x) and x > 0 for x in (a.days, a.dt, a.save_hours)):
        ap.error('days, dt, save-hours must be finite and > 0')
    total = a.days * b3.SEC_DAY / a.dt
    save = a.save_hours * 3600. / a.dt
    if abs(total - round(total)) > 1e-8 or abs(save - round(save)) > 1e-8:
        ap.error('days*86400/dt and save-hours*3600/dt must be integers')
    a.steps, a.save_every = int(round(total)), int(round(save))
    try:
        a.lights = [float(x) for x in a.lights.split(',')]
    except ValueError:
        ap.error('--lights must be comma-separated numbers')
    if not a.lights or not all(math.isfinite(x) and x >= 0 for x in a.lights):
        ap.error('invalid --lights')
    a.treatments = (('replete',) + SELECTED_TREATMENTS
                    if a.include_replete else SELECTED_TREATMENTS)
    return a


def parameter(model, name):
    found = [p for p in model.parameters if p.name == 'P1/' + name]
    if len(found) != 1:
        raise RuntimeError(f'Expected exactly one P1/{name} parameter, got {len(found)}')
    return float(found[0].value)


def quota_limitation(n_to_c, p_to_c, params):
    qnlc, qplc, xqcn, xqcp = params
    n_denom = xqcn * QN_REDFIELD - qnlc
    p_denom = xqcp * QP_REDFIELD - qplc
    if min(n_denom, p_denom) <= 0:
        raise ValueError('Invalid positive quota-limitation denominators')
    iN = float(np.clip((n_to_c - qnlc) / n_denom, 0, 1))
    iP = float(np.clip((p_to_c - qplc) / p_denom, 0, 1))
    if iN > 1-1e-9 and iP > 1-1e-9:
        active = 'both_saturated'
    elif abs(iN - iP) < 1e-7:
        active = 'tie'
    else:
        active = 'N' if iN < iP else 'P'
    return iN, iP, active


def check_saved(row, model, rates, case):
    if not np.isfinite(rates).all():
        raise RuntimeError(f'Nonfinite rates in {case}')
    for col in ('P1_gross_C_fixation', 'P1_NO3_uptake',
                'P1_NH4_uptake', 'P1_PO4_uptake'):
        if not math.isfinite(float(row[col])):
            raise RuntimeError(f'Missing or invalid flux diagnostic {col} in {case}')
    if row['condJ'] < 0 or not math.isfinite(row['condJ']):
        raise RuntimeError(f'Invalid condJ in {case}')


def integrate(a, sw, tr, mode):
    m, initial, _, _, fdep, swdep, parvar = b3.build_model()
    nidx = {v.name: i for i, v in enumerate(m.state_variables)}
    for name in b3.STATE_NAMES:
        if name not in nidx:
            raise RuntimeError(f'Missing state {name}')
    for name in SOURCE_VERIFIED_SIGNED:
        if name not in nidx:
            raise RuntimeError(f'Missing source-verified signed state {name}')
    params = tuple(parameter(m, x) for x in ('qnlc', 'qplc', 'xqcn', 'xqcp'))
    if int(parameter(m, 'Limnut')) != 1:
        raise RuntimeError('This branch reconstruction assumes P1/Limnut=1 (hard minimum)')
    bounds = (parameter(m, 'fC_min'), parameter(m, 'fC_max'))
    rate_cap = parameter(m, 'fC_rate_max')
    if not np.allclose(bounds, (b3.FC_MIN, b3.FC_MAX), rtol=0, atol=1e-12):
        raise RuntimeError(f'Allocation bounds changed: {bounds}')
    if abs(rate_cap - 1.5) > 1e-10:
        raise RuntimeError(f'Unexpected fC rate cap: {rate_cap}')

    diag = b3.diagnose(m, b3.REQUIRED_DIAG + b3.OPTIONAL_DIAG)
    for d in b3.REQUIRED_DIAG + NEEDED_FLUXES:
        if diag[d] is None:
            raise RuntimeError(f'Missing required P1 diagnostic: {d}')
    bathN, bathP = b3.TREATMENTS[tr]
    bath = {'N3/n': bathN, 'N4/n': b3.NH4_BATH,
            'N1/p': bathP, 'N5/s': b3.SI_BATH}
    bi = np.array([nidx[x] for x in bath], dtype=int)
    bv = np.array(list(bath.values()), dtype=float)
    state = initial.copy()
    state[bi] = bv
    fc = b3.F0
    swdep.value = float(sw)
    # The Python strict nonnegativity check retains the same single
    # source-verified exception as the successful B3_v2 experiment.
    signed_idx = {nidx[s] for s in SOURCE_VERIFIED_SIGNED}
    nonneg_idx = np.array([i for i in np.flatnonzero(state >= 0)
                           if i not in signed_idx], dtype=int)
    signed_min = float(state[nidx['G2/o_deep']])
    history = []
    flux_names = ('P1_gross_C_fixation', 'P1_NO3_uptake',
                  'P1_NH4_uptake', 'P1_PO4_uptake')
    integrated = {x: 0.0 for x in flux_names}
    mean_specific = {x: 0.0 for x in flux_names}
    active_counts = {x: 0 for x in ('N', 'P', 'tie', 'both_saturated')}
    switches = 0
    last_active = None
    bounds_hits = 0
    minimum_cond = float('inf')
    for k in range(a.steps + 1):
        row, rates = b3.sample(m, state, fc, sw, k, a.dt, diag, parvar,
                              tr, bath, nidx, fdep, b3.PAR_UMOL_PER_J)
        case = f'{tr} SW={sw} {mode} step={k}'
        check_saved(row, m, rates, case)
        iN, iP, active = quota_limitation(row['N_to_C_native'],
                                           row['P_to_C_native'], params)
        row.update(mode=mode, i_N=iN, i_P=iP, i_NP=min(iN, iP),
                   limiting_branch=active)
        minimum_cond = min(minimum_cond, row['condJ'])
        if k % a.save_every == 0 or k == a.steps:
            history.append(row)
        if k == a.steps:
            break
        active_counts[active] += 1
        if last_active is not None and active != last_active:
            switches += 1
        last_active = active
        dt_day = a.dt / b3.SEC_DAY
        for metric in flux_names:
            integrated[metric] += dt_day * row[metric]
            mean_specific[metric] += dt_day * row[metric] / row['P1_C']
        next_state = state + a.dt * rates
        next_state[bi] = bv  # explicit open reservoir
        signed_min = min(signed_min, next_state[nidx['G2/o_deep']])
        if not np.isfinite(next_state).all():
            raise RuntimeError(f'Nonfinite state: {case}')
        negatives = nonneg_idx[next_state[nonneg_idx] < 0]
        if len(negatives):
            j = int(negatives[np.argmin(next_state[negatives])])
            raise RuntimeError(f'Unexpected negative state {m.state_variables[j].name}='
                               f'{next_state[j]:.6g} in {case}; do not clip')
        if mode == 'dynamic':
            trial = fc + dt_day * row['dfCdt_d_inv']
            new_fc = float(np.clip(trial, *bounds))
            bounds_hits += int(trial != new_fc)
        elif mode == 'fixed_ref':
            # DO NOT turn off the machinery multipliers; freeze only host fC.
            new_fc = b3.F0
        else:
            raise ValueError(f'Unknown mode: {mode}')
        state, fc = next_state, new_fc

    first, last = history[0], history[-1]
    tail = [r for r in history if r['t_days'] >= .75*a.days]
    summary = dict(treatment=tr, shortwave_W_m2=sw, mode=mode,
                   duration_days=a.days, dt_seconds=a.dt,
                   initial_P1_C=first['P1_C'], final_P1_C=last['P1_C'],
                   final_fC=last['fC'], fC_drift_last_quarter=(
                       abs(tail[-1]['fC']-tail[0]['fC']) if len(tail)>1 else float('nan')),
                   final_N_to_C=last['N_to_C_native'],
                   final_P_to_C=last['P_to_C_native'],
                   final_i_N=last['i_N'], final_i_P=last['i_P'],
                   final_limiter=last['limiting_branch'],
                   percent_time_N=100*active_counts['N']/a.steps,
                   percent_time_P=100*active_counts['P']/a.steps,
                   percent_time_saturated=100*active_counts['both_saturated']/a.steps,
                   branch_switch_count=switches,
                   fC_bound_clamp_count=bounds_hits,
                   min_condJ=minimum_cond,
                   min_G2_o_deep_signed=signed_min)
    for metric in flux_names:
        summary['integral_days_times_' + metric] = integrated[metric]
        summary['mean_specific_' + metric] = mean_specific[metric] / a.days
    summary['mean_specific_total_inorganic_N_uptake'] = (
        mean_specific['P1_NO3_uptake'] + mean_specific['P1_NH4_uptake']) / a.days
    return history, summary


def dump_csv(path, rows):
    if path.exists():
        raise FileExistsError(path)
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def delta_pct(dynamic, fixed):
    if not np.isfinite(dynamic) or not np.isfinite(fixed) or abs(fixed) < 1e-14:
        return float('nan')
    return 100 * (dynamic-fixed)/abs(fixed)


def plot_light(out, sw, case_rows, treatments):
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.6), constrained_layout=True)
    definitions = [
        ('fC', r'Carbon allocation $f_C$'),
        ('gross_C_fixation_per_C', 'Gross fixation / P1 carbon (native diagnostic)'),
        ('P1_NO3_uptake', 'NO3 uptake flux (native diagnostic)'),
        ('P1_PO4_uptake', 'PO4 uptake flux (native diagnostic)'),
    ]
    for ax, (key, ylabel) in zip(axes.flat, definitions):
        for tr in treatments:
            color = b3.C_COLORS[tr]
            for mode, style in [('dynamic', '-'), ('fixed_ref', '--')]:
                r = case_rows[(tr, mode)]
                y = np.array([x[key] for x in r], dtype=float)
                if key in ('P1_NO3_uptake', 'P1_PO4_uptake'):
                    y = y / np.array([x['P1_C'] for x in r], dtype=float)
                ax.plot([x['t_days'] for x in r], y, ls=style,
                        color=color, lw=1.6 if mode == 'fixed_ref' else 2,
                        label=f'{tr.replace("_", " ")} {"dynamic" if mode=="dynamic" else "fixed"}')
        ax.set(xlabel='Elapsed days', ylabel=(ylabel + (' / P1 carbon' if key in (
            'P1_NO3_uptake', 'P1_PO4_uptake') else '')))
        ax.grid(alpha=0.18)
        if key == 'fC':
            ax.axhline(b3.F0, lw=.8, ls=':', color='0.5')
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc='lower center',
               bbox_to_anchor=(0.5, -0.04), ncol=3, fontsize=9)
    fig.suptitle(f'P1 physiological response: dynamic versus reference-fixed allocation; '
                 f'SW={sw:g} W m$^{{-2}}$', fontsize=13)
    for ext in ('png', 'pdf'):
        p = out / f'B4_paired_SW{sw:g}.{ext}'
        if p.exists():
            raise FileExistsError(p)
        fig.savefig(p, dpi=240, bbox_inches='tight')
    plt.close(fig)


def main():
    a = parse_args()
    if a.output.exists() and any(a.output.iterdir()):
        raise FileExistsError(f'Choose a new output directory: {a.output}')
    a.output.mkdir(parents=True, exist_ok=True)
    print('ERSEM-ALLOC B4: P1 dynamic vs fixed-reference allocation', flush=True)
    print('Time step', a.dt, 'seconds; duration', a.days, 'days;',
          'lights', a.lights, 'treatments', a.treatments, flush=True)
    print('External NO3, NH4, PO4, Si are fixed baths (NOT closed budgets).')
    print('Fixed-reference arm retains all fC capacity multipliers; only fC memory is frozen.')
    print('Replete excluded by default due to B3 gradient switching. No originals modified.')
    pair_summaries = []
    all_summaries = []
    for sw in a.lights:
        series = {}
        for tr in a.treatments:
            arms = {}
            for mode in ('fixed_ref', 'dynamic'):
                rows, s = integrate(a, sw, tr, mode)
                swlabel = f'{sw:g}'.replace('.', 'p')
                dump_csv(a.output / f'B4_SW{swlabel}_{tr}_{mode}.csv', rows)
                series[(tr, mode)] = rows
                arms[mode] = s
                all_summaries.append(s)
                print(f'{sw:g} W/m2 {tr:8s} {mode:10s}: '
                      f'fC={s["final_fC"]:.6f}, P1C={s["final_P1_C"]:.5g}, '
                      f'iN={s["final_i_N"]:.4f}, iP={s["final_i_P"]:.4f}, '
                      f'limiter={s["final_limiter"]}, branch-switches={s["branch_switch_count"]}',
                      flush=True)
            # Same initial state and environmental settings across arms.
            a0 = series[(tr, 'dynamic')][0]
            b0 = series[(tr, 'fixed_ref')][0]
            for key in ('P1_C', 'P1_N', 'P1_P', 'P1_Chl', 'fC',
                        'P1_gross_C_fixation', 'P1_NO3_uptake', 'P1_PO4_uptake'):
                if not np.isclose(a0[key], b0[key], rtol=0, atol=1e-11):
                    raise RuntimeError(f'Paired initial-state mismatch {tr} {sw}: {key}')
            dyn, fix = arms['dynamic'], arms['fixed_ref']
            pair = {'shortwave_W_m2': sw, 'treatment': tr, 'days': a.days,
                    'dt_seconds': a.dt, 'fC_dynamic_final': dyn['final_fC'],
                    'fC_fixed_final': fix['final_fC'],
                    'limiter_dynamic_final': dyn['final_limiter'],
                    'limiter_fixed_final': fix['final_limiter'],
                    'P1_C_dynamic_final': dyn['final_P1_C'],
                    'P1_C_fixed_final': fix['final_P1_C'],
                    'P1_C_percent_change_dynamic_vs_fixed': delta_pct(
                        dyn['final_P1_C'], fix['final_P1_C'])}
            for metric in ('P1_gross_C_fixation', 'P1_NO3_uptake',
                           'P1_NH4_uptake', 'P1_PO4_uptake'):
                k = 'mean_specific_' + metric
                pair[metric + '_specific_dynamic'] = dyn[k]
                pair[metric + '_specific_fixed'] = fix[k]
                pair[metric + '_specific_pct_change'] = delta_pct(dyn[k], fix[k])
            nk = 'mean_specific_total_inorganic_N_uptake'
            pair['total_inorganic_N_specific_dynamic'] = dyn[nk]
            pair['total_inorganic_N_specific_fixed'] = fix[nk]
            pair['total_inorganic_N_specific_pct_change'] = delta_pct(dyn[nk], fix[nk])
            pair_summaries.append(pair)
            print(f' PAIRED {sw:g} {tr}: specific carbon fixation '
                  f'{pair["P1_gross_C_fixation_specific_pct_change"]:+.2f}%, '
                  f'inorg N uptake {pair["total_inorganic_N_specific_pct_change"]:+.2f}%, '
                  f'PO4 uptake {pair["P1_PO4_uptake_specific_pct_change"]:+.2f}%, '
                  f'final P1 carbon {pair["P1_C_percent_change_dynamic_vs_fixed"]:+.2f}%',
                  flush=True)
        plot_light(a.output, sw, series, a.treatments)
    dump_csv(a.output / 'B4_all_arm_summary.csv', all_summaries)
    dump_csv(a.output / 'B4_paired_comparison_summary.csv', pair_summaries)
    meta = dict(experiment='B4 dynamic-vs-fixed P1 reference allocation',
                duration_days=a.days, dt_seconds=a.dt, lights=a.lights,
                treatments=a.treatments, reference_fC=b3.F0,
                free_internal_states=True,
                bath={'NO3': 'specified by treatment', 'NH4': b3.NH4_BATH,
                      'PO4': 'specified by treatment', 'Si': b3.SI_BATH},
                no_original_input_changes=True, numerics='forward Euler; not NEMO',
                signed_benthic_state='G2/o_deep allowed negative by source')
    with (a.output / 'B4_experiment_metadata.json').open('w') as f:
        json.dump(meta, f, indent=2)
    print('Summary:', a.output / 'B4_paired_comparison_summary.csv', flush=True)
    print('Do not infer optimal allocation, fitness advantage, or natural NP limitation')
    print('without checking the branch diagnostics, flux definitions and timestep convergence.')


if __name__ == '__main__':
    main()
