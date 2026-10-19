#!/usr/bin/env python3
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
"""Two-panel supplementary physiology figure, from regenerated CSVs only.

(a) Instantaneous Chl:C specific-rate balance at four prescribed fC values and
    digitised Geider et al. 1997 Fig 3A data (observational context).
(b) Gross carbon fixation per P1 carbon, fixed Chl:C across treatments,
    normalised to the fC=0.24 model response at PAR=300 umol photons m-2 s-1.
No time integration, fitting, dynamic-fC optimisation, or validation claim.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/'results';OBS=ROOT/'observations/geider1997_thalassiosira_observations.csv'
FIGURES=ROOT/'figures'
FC=(0.12,0.24,0.38,0.50)
COLORS={0.12:'#327B97',0.24:'#2255A4',0.38:'#CF7833',0.50:'#8B4F90'}

def read(path):
    if not path.exists():raise FileNotFoundError(path)
    with path.open(newline='') as f:return list(csv.DictReader(f))

def cols(rows,x,y):
    return np.asarray([float(r[x]) for r in rows]),np.asarray([float(r[y]) for r in rows])

def main():
    obs=read(OBS);bal=read(RESULTS/'ersem_alloc_fixed_fc_sensitivity.csv')
    pi=read(RESULTS/'ersem_alloc_p1_photosynthesis_irradiance.csv')
    if len(obs)!=11:raise ValueError('Expected 11 digitised Geider observations')
    if any(not r['status']=='ONE_ROOT' for r in bal):
        raise ValueError('Missing or multiple Chl:C roots; review gaps before plotting')
    if not pi:raise ValueError('No P-I results')
    xo,yo=cols(obs,'PAR_umol_photons_m2_s','ChlC_g_g')
    plt.rcParams.update({'font.family':'sans-serif','font.size':9.3,
        'axes.labelsize':10,'xtick.labelsize':9,'ytick.labelsize':9,
        'legend.fontsize':8.0,'pdf.fonttype':42,'ps.fonttype':42,
        'axes.linewidth':0.85})
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(9.3,3.85),layout='constrained')
    for ax in (ax1,ax2):
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
        ax.tick_params(direction='out',length=3.8,width=0.8)
        ax.set_xlim(0,600)
    # Panel a: all roots 48/48 per fC; same model state other than P1 Chl and fC.
    for fc in FC:
        rr=[r for r in bal if abs(float(r['fC'])-fc)<1e-9]
        if len(rr)!=48:raise ValueError(f'Expected 48 points for fC={fc}; found {len(rr)}')
        x,y=cols(rr,'model_PAR_umol_photons_m2_s','balanced_ChlC_g_g')
        order=np.argsort(x)
        ax1.plot(x[order],y[order],color=COLORS[fc],lw=2.3 if fc==0.24 else 1.8,
                 linestyle='-' if fc==0.24 else '--',label=rf'$f_C={fc:.2f}$')
    ax1.scatter(xo,yo,s=25,marker='o',color='black',zorder=8,label='Geider et al. (1997)')
    ax1.set_ylim(0,0.066)
    ax1.set_xlabel(r'PAR ($\mu$mol photons m$^{-2}$ s$^{-1}$)')
    ax1.set_ylabel(r'Balanced Chl:$C$ (g Chl $a$ g$^{-1}$ C)')
    ax1.legend(frameon=False,loc='upper right',handlelength=2.1)
    # Panel b: isolate fC by holding Chl:C = 0.025, all other states unchanged.
    pis=[r for r in pi if r['mode']=='FIXED_CHLC']
    reference=[r for r in pis if abs(float(r['fC'])-0.24)<1e-9]
    xr,yr=cols(reference,'PAR_umol_photons_m2_s','gross_C_fixation_per_P1C_native')
    order=np.argsort(xr)
    scale=float(np.interp(300.,xr[order],yr[order]))
    if not np.isfinite(scale) or scale<=0:raise ValueError('Invalid fC=0.24 / PAR300 normalization')
    for fc in FC:
        rr=[r for r in pis if abs(float(r['fC'])-fc)<1e-9]
        if len(rr)<20:raise ValueError(f'Insufficient P-I points for fC={fc}')
        x,y=cols(rr,'PAR_umol_photons_m2_s','gross_C_fixation_per_P1C_native')
        order=np.argsort(x)
        ax2.plot(x[order],y[order]/scale,color=COLORS[fc],lw=2.3 if fc==0.24 else 1.8,
                 linestyle='-' if fc==0.24 else '--',label=rf'$f_C={fc:.2f}$')
    ax2.axhline(1.0,lw=0.65,linestyle=':',color='#888888')
    ax2.set_ylim(bottom=0)
    ax2.set_xlabel(r'PAR ($\mu$mol photons m$^{-2}$ s$^{-1}$)')
    ax2.set_ylabel('P1 gross C fixation / biomass\n(normalised to reference at PAR = 300)')
    ax2.legend(frameon=False,loc='upper left',handlelength=2.1)
    ax1.text(-0.095,1.025,'a',fontweight='bold',fontsize=12,transform=ax1.transAxes)
    ax2.text(-0.095,1.025,'b',fontweight='bold',fontsize=12,transform=ax2.transAxes)
    FIGURES.mkdir(exist_ok=True)
    for ext in ('png','pdf','svg'):
        dest=FIGURES/f'supp_geider_pi.{ext}'
        fig.savefig(dest,dpi=400,bbox_inches='tight',facecolor='white')
        print('Saved:',dest)
    plt.close(fig)
    print('Note: the normalised P-I values are NOT measured photosynthesis, '
          'and the Chl:C model curves are conditional roots, NOT time-integrated fC dynamics.')

if __name__=='__main__': main()
