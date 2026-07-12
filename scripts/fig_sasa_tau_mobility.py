#!/usr/bin/env python
"""Paper figure: T_k and sigma_k vs SASA, colored by dipole mobility. 2×3 layout.

For appendix. Uses corrected tri-exp fit (matching Table 1). Top row: relaxation
timescales T_k; bottom row: per-component amplitudes sigma_k = sigma_m * sqrt(f_k).
"""
from __future__ import annotations
import csv, json
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

mpl.use('Agg')
plt.rcParams.update({
    'figure.dpi': 130, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'font.size': 13, 'axes.labelsize': 14, 'axes.titlesize': 13,
    'xtick.labelsize': 11, 'ytick.labelsize': 11,
    'axes.grid': True, 'grid.alpha': 0.25, 'axes.axisbelow': True,
    'font.family': 'DejaVu Sans',
})

RESULTS = Path(__file__).resolve().parent.parent / 'results'
OUT = RESULTS / 'structural_correlation'

SASA = np.array([0.194, 0.703, 0.046, 0.488, 0.029, 0.066, 0.157, 1.003])
TRP_NAMES = [f'Trp{i+1}' for i in range(8)]
COL_TITLES = ['Solvent libration', 'Water reorientation', 'Protein conformational']


def load_fit_params():
    """Load (tau_k, f_k) from corrected_triexp.csv. Returns taus (8,3), fks (8,3)."""
    taus = np.zeros((8, 3))
    fks = np.zeros((8, 3))
    with open(RESULTS / 'paper_figures' / 'corrected_triexp.csv') as f:
        for row in csv.DictReader(f):
            i = int(row['site'].replace('Trp', '')) - 1
            taus[i] = [float(row['tau1_ps']), float(row['tau2_ps']),
                       float(row['tau3_ps'])]
            fks[i] = [float(row['A1']), float(row['A2']), float(row['A3'])]
    return taus, fks


def load_sigma():
    """Load per-site total sigma (cm^-1) from table1.csv (skips the mean row)."""
    sigma = np.zeros(8)
    with open(RESULTS / 'paper_figures' / 'table1.csv') as f:
        for row in csv.DictReader(f):
            if row['site'] == 'mean':
                continue
            i = int(row['site'].replace('Trp', '')) - 1
            sigma[i] = float(row['sigma_total_cm'])
    return sigma


def load_mobility():
    with open(OUT / 'indole_mobility.json') as f:
        return np.array(json.load(f)['angular_deviation_local_deg'])


taus, fks = load_fit_params()
sigma_total = load_sigma()
theta = load_mobility()

# per-component amplitude: sigma_{m,k} = sigma_m * sqrt(f_k)
sigma_k = sigma_total[:, None] * np.sqrt(fks)   # (8, 3)

T1_fs = taus[:, 0] * 1000
T2_ps = taus[:, 1]
T3_ns = taus[:, 2] / 1000

row1 = [(T1_fs, r'$T_1$ (fs)'),
        (T2_ps, r'$T_2$ (ps)'),
        (T3_ns, r'$T_3$ (ns)')]
row2 = [(sigma_k[:, 0], r'$\sigma_1$ (cm$^{-1}$)'),
        (sigma_k[:, 1], r'$\sigma_2$ (cm$^{-1}$)'),
        (sigma_k[:, 2], r'$\sigma_3$ (cm$^{-1}$)')]

fig, axs = plt.subplots(2, 3, figsize=(12, 7.4), sharex=True,
                        gridspec_kw={'hspace': 0.14, 'wspace': 0.30})

norm = mpl.colors.Normalize(vmin=theta.min(), vmax=theta.max())
cmap = plt.cm.RdYlBu_r

for col in range(3):
    for row, (vals, ylabel) in enumerate([row1[col], row2[col]]):
        ax = axs[row, col]
        sc = ax.scatter(SASA, vals, c=theta, cmap=cmap, norm=norm,
                        s=80, zorder=5, edgecolors='black', linewidths=0.5)
        for i in range(8):
            ax.annotate(TRP_NAMES[i], (SASA[i], vals[i]),
                        textcoords='offset points', xytext=(6, 3), fontsize=7.5)
        ax.set_ylabel(ylabel)
        ax.set_xlim(-0.12, 1.17)
        if row == 0:
            ax.set_title(COL_TITLES[col], fontsize=11)
        if row == 1:
            ax.set_xlabel('SASA (nm$^2$)', fontsize=12)

# hide x tick labels on the top row (shared axis)
for col in range(3):
    plt.setp(axs[0, col].get_xticklabels(), visible=False)

# row 2 (sigma): unify y-range across the three panels
ylo = min(axs[1, c].get_ylim()[0] for c in range(3))
yhi = max(axs[1, c].get_ylim()[1] for c in range(3))
for c in range(3):
    axs[1, c].set_ylim(ylo-50, yhi+50)
    axs[1, c].set_yticks(np.arange(0, yhi+1, 200))
    axs[1, c].set_yticks(np.arange(0, yhi+1, 100), minor=True)
    fmt = ScalarFormatter(useMathText=True)
    fmt.set_scientific(True)
    fmt.set_powerlimits((0, 0))
    axs[1, c].yaxis.set_major_formatter(fmt)
axs[1, 0].yaxis.get_offset_text().set_fontsize(10)

cbar = fig.colorbar(sc, ax=axs, orientation='horizontal',
                    fraction=0.025, pad=0.12, aspect=36)
cbar.set_label(r'$\theta_{\mathrm{loc}}$ (deg)', fontsize=12)

for ext in ['png', 'pdf']:
    fig.savefig(OUT / f'fig_sasa_tau_mobility.{ext}')
plt.close(fig)
print(f"Saved -> {OUT}/fig_sasa_tau_mobility.{{png,pdf}}")
