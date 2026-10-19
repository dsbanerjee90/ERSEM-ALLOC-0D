# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
#
# Section 2 / B4 figure for the manuscript.
#
# Purpose
# -------
# Build a clean four-panel manuscript figure from the B4 paired dynamic-vs-fixed
# simulations, using DAILY-MEAN trajectories for the time-series panels so that
# short-timescale numerical transients do not dominate the visual message.
#
# Inputs expected in --input directory:
#   - B4_paired_comparison_summary.csv
#   - B4_SW20_low_N_dynamic.csv
#   - B4_SW20_low_N_fixed_ref.csv
#   - B4_SW20_low_P_dynamic.csv
#   - B4_SW20_low_P_fixed_ref.csv
#   - B4_SW100_low_N_dynamic.csv
#   - B4_SW100_low_N_fixed_ref.csv
#   - B4_SW100_low_P_dynamic.csv
#   - B4_SW100_low_P_fixed_ref.csv
#
# Outputs written into --output directory:
#   - B4_scientific_main_dailymean.png
#   - B4_scientific_main_dailymean.pdf
#   - B4_scientific_main_dailymean.svg
#
# Example:
# python plot_B4_main_dailymean.py \
#   --input ../results/B4_7d_dt150 \
#   --output ../figures/B4_main_dailymean

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

TREATMENTS = ('low_N', 'low_P')
LIGHTS = (20.0, 100.0)
EXPECTED_LIMITER = {'low_N': 'N', 'low_P': 'P'}
SUMMARY_REQUIRED = (
    'shortwave_W_m2', 'treatment',
    'fC_dynamic_final', 'fC_fixed_final',
    'limiter_dynamic_final',
    'P1_gross_C_fixation_specific_pct_change',
    'total_inorganic_N_specific_pct_change',
    'P1_PO4_uptake_specific_pct_change',
    'P1_C_percent_change_dynamic_vs_fixed',
    'P1_C_dynamic_final', 'P1_C_fixed_final',
)


def read_csv(path: Path):
    with path.open(newline='') as fp:
        rd = csv.DictReader(fp)
        names = list(rd.fieldnames or [])
        if not names:
            raise ValueError(f'No CSV columns in {path}')
        rows = list(rd)
    if not rows:
        raise ValueError(f'Empty CSV: {path}')
    return names, rows


def read_summary(path: Path):
    names, rows = read_csv(path)
    missing = set(SUMMARY_REQUIRED) - set(names)
    if missing:
        raise ValueError(f'Missing required summary columns: {sorted(missing)}')

    out = {}
    for r in rows:
        tr = r['treatment']
        sw = float(r['shortwave_W_m2'])
        if tr not in TREATMENTS or sw not in LIGHTS:
            continue
        key = (sw, tr)
        out[key] = r
        if r['limiter_dynamic_final'] != EXPECTED_LIMITER[tr]:
            raise ValueError(f'{key} has unexpected limiter {r["limiter_dynamic_final"]}')
        if abs(float(r['fC_fixed_final']) - 0.24) > 1e-8:
            raise ValueError(f'Fixed reference is not fC=0.24 for {key}')

    for sw in LIGHTS:
        for tr in TREATMENTS:
            if (sw, tr) not in out:
                raise ValueError(f'Missing summary row for {(sw, tr)}')
    return out


def read_trajectories(folder: Path, summary_rows):
    required = ('t_days', 'fC', 'P1_C')
    curves = {}
    for sw in LIGHTS:
        for tr in TREATMENTS:
            for mode in ('dynamic', 'fixed_ref'):
                path = folder / f'B4_SW{sw:g}_{tr}_{mode}.csv'
                if not path.is_file():
                    raise FileNotFoundError(path)
                names, rows = read_csv(path)
                missing = set(required) - set(names)
                if missing:
                    raise ValueError(f'{path} missing {sorted(missing)}')
                t = np.asarray([float(r['t_days']) for r in rows], dtype=float)
                fc = np.asarray([float(r['fC']) for r in rows], dtype=float)
                pc = np.asarray([float(r['P1_C']) for r in rows], dtype=float)
                if abs(t[0]) > 1e-8 or abs(t[-1] - 7.0) > 1e-6:
                    raise ValueError(f'Unexpected time span in {path}')
                if mode == 'fixed_ref' and np.max(np.abs(fc - 0.24)) > 1e-7:
                    raise ValueError(f'Fixed arm is not constant in {path}')
                expected_fc = float(summary_rows[(sw, tr)]['fC_dynamic_final' if mode == 'dynamic' else 'fC_fixed_final'])
                if abs(fc[-1] - expected_fc) > 1e-6:
                    raise ValueError(f'fC mismatch between trajectory and summary in {path}')
                curves[(sw, tr, mode)] = (t, fc, pc)
    return curves


def compute_daily_means(t, y, ndays=7):
    """Daily means of seven non-overlapping 24 h bins [0,1), ..., [6,7).

    Bin 1 is displayed at x=1 (experiment day 1), ..., bin 7 at x=7.
    The exact endpoint t=7 is not included in a daily mean; it is used
    separately for the final-state biomass comparison in panel d.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    x_daily = []
    y_daily = []
    for day in range(ndays):
        mask = (t >= day) & (t < day + 1)
        if not np.any(mask):
            raise ValueError(f'No samples found for day bin {day}')
        x_daily.append(day + 1.0)
        y_daily.append(float(np.mean(y[mask])))
    return np.asarray(x_daily), np.asarray(y_daily)


def get_colors():
    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    return {'low_N': cycle[0], 'low_P': cycle[1]}


def make_axes_pretty(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(False)
    ax.set_axisbelow(True)


def mark_panel(ax, letter, title):
    # The descriptive panel titles belong in the manuscript caption.
    ax.text(-0.10, 1.025, letter, transform=ax.transAxes,
            fontsize=9.5, fontweight='bold', ha='left', va='bottom')


def limiting_uptake_pct(row, treatment):
    if treatment == 'low_N':
        return float(row['total_inorganic_N_specific_pct_change'])
    return float(row['P1_PO4_uptake_specific_pct_change'])


def save_all(fig, outdir: Path, stem: str):
    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ('png', 'pdf', 'svg'):
        path = outdir / f'{stem}.{ext}'
        fig.savefig(path, dpi=450, bbox_inches='tight', facecolor='white')
        print(f'Saved: {path}')
    plt.close(fig)


def plot_main_figure(outdir: Path, summary_rows, curves):
    colors = get_colors()
    fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.65), layout='constrained')

    # Panel a and b: DAILY-MEAN dynamic fC trajectories at low and high light.
    panel_specs = [('a', 20.0, 'Lower light (20 W m$^{-2}$)'),
                   ('b', 100.0, 'Higher light (100 W m$^{-2}$)')]

    for col, (letter, sw, title) in enumerate(panel_specs):
        a = ax[0, col]
        make_axes_pretty(a)
        mark_panel(a, letter, title)
        # Experimental conditions, not descriptive subplot titles.
        a.text(0.5, 1.015, rf'{int(sw)} W m$^{{-2}}$',
               transform=a.transAxes, fontsize=8.2,
               ha='center', va='bottom')

        for tr, text in [('low_N', 'N-limited'), ('low_P', 'P-limited')]:
            t, fc, _ = curves[(sw, tr, 'dynamic')]
            xd, yd = compute_daily_means(t, fc)
            a.plot(xd, yd, lw=1.8, color=colors[tr], label=text)

        # Fixed reference line remains constant; show as horizontal reference.
        a.axhline(0.24, lw=1.2, linestyle='--', color='#666666', label=r'Fixed reference ($f_C = 0.24$)')

        a.set_xlim(0.85, 7.15)
        a.set_xticks(np.arange(1, 8, 1))
        a.set_ylim(0.20, 0.55)
        a.set_xlabel('Experiment day')
        if col == 0:
            a.set_ylabel(r'Effective carbon allocation, $f_C$')

    handles, labels = ax[0, 0].get_legend_handles_labels()
    ax[0, 0].legend(handles, labels, loc='upper left', frameon=False)

    # Panel c: phase plane of fixation gain vs limiting nutrient uptake change.
    c = ax[1, 0]
    make_axes_pretty(c)
    mark_panel(c, 'c', 'Photosynthesis–nutrient uptake trade-off')
    c.axvline(0, linestyle='--', lw=0.8, alpha=0.65, color='#777777')
    c.axhline(0, linestyle='--', lw=0.8, alpha=0.65, color='#777777')

    for tr, short in [('low_N', 'N'), ('low_P', 'P')]:
        for sw in LIGHTS:
            r = summary_rows[(sw, tr)]
            x = float(r['P1_gross_C_fixation_specific_pct_change'])
            y = limiting_uptake_pct(r, tr)
            c.scatter(x, y, s=38, color=colors[tr], zorder=4)
            c.annotate(f'{short} · {int(sw)}', (x, y), xytext=(5, 5),
                       textcoords='offset points', fontsize=7.2)

    c.set_xlim(-3, 52)
    c.set_ylim(-73, 5)
    c.set_xlabel(r'$\Delta$ gross C fixation / biomass (%)')
    c.set_ylabel(r'$\Delta$ limiting nutrient uptake / biomass (%)')

    # Panel d: final carbon comparison.
    d = ax[1, 1]
    make_axes_pretty(d)
    mark_panel(d, 'd', 'Less phytoplankton carbon at day seven')
    d.grid(False)

    rowspec = [(20.0, 'low_N'), (20.0, 'low_P'), (100.0, 'low_N'), (100.0, 'low_P')]
    yvals = [3, 2, 1, 0]
    ylabels = [
        'N-limited · 20 W m$^{-2}$',
        'P-limited · 20 W m$^{-2}$',
        'N-limited · 100 W m$^{-2}$',
        'P-limited · 100 W m$^{-2}$',
    ]

    xmax = 0.0
    for y, (sw, tr) in zip(yvals, rowspec):
        r = summary_rows[(sw, tr)]
        dynamic = float(r['P1_C_dynamic_final'])
        fixed = float(r['P1_C_fixed_final'])
        delta = float(r['P1_C_percent_change_dynamic_vs_fixed'])
        xmax = max(xmax, dynamic, fixed)
        d.plot([dynamic, fixed], [y, y], lw=1.5, alpha=0.75, color=colors[tr])
        d.scatter([fixed], [y], marker='D', s=36, facecolors='none', edgecolors=colors[tr], linewidths=1.1, zorder=4)
        d.scatter([dynamic], [y], marker='o', s=46, color=colors[tr], zorder=5)
        d.annotate(f'{delta:+.1f}%', (max(dynamic, fixed), y), xytext=(5, 0),
                   textcoords='offset points', va='center', fontsize=7.1)

    d.set_yticks(yvals)
    d.set_yticklabels(ylabels)
    d.set_xlim(0, xmax * 1.34)
    d.set_ylim(-0.3, 3.95)  # Reserve clear space for the two-line legend.
    d.set_xlabel('P1 carbon at day 7 (mg C m$^{-3}$)')

    legend_handles = [
        Line2D([], [], marker='D', linestyle='none', markerfacecolor='none', markeredgecolor='black', markersize=5.0, label='Fixed reference'),
        Line2D([], [], marker='o', linestyle='none', markerfacecolor='black', markeredgecolor='black', markersize=5.0, label='Dynamic allocation'),
    ]
    d.legend(handles=legend_handles, frameon=False, loc='upper right')

    # No figure-level title: report the mechanism in the paper caption.
    save_all(fig, outdir, 'fig3')


def main():
    p = argparse.ArgumentParser(description='Create the daily-mean B4 scientific figure.')
    p.add_argument('--input', required=True, type=Path,
                   help='B4 result directory containing B4_paired_comparison_summary.csv and trajectory CSVs')
    p.add_argument('--output', required=True, type=Path,
                   help='Output directory for figure files')
    args = p.parse_args()

    folder = args.input
    if not folder.is_dir():
        raise FileNotFoundError(folder)

    summary_path = folder / 'B4_paired_comparison_summary.csv'
    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)

    args.output.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.size': 7.6,
        'axes.labelsize': 8.3,
        'axes.titlesize': 8.3,
        'xtick.labelsize': 7.1,
        'ytick.labelsize': 7.1,
        'legend.fontsize': 7.2,
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.major.size': 3.5,
        'ytick.major.size': 3.5,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    summary_rows = read_summary(summary_path)
    curves = read_trajectories(folder, summary_rows)
    plot_main_figure(args.output, summary_rows, curves)


if __name__ == '__main__':
    main()