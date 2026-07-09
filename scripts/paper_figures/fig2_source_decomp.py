#!/usr/bin/env python
"""Fig 2: source decomposition + dielectric screening.

Left:  grouped sigma bars (protein/water/nucleotide/ions) per site, overlaid
        with total sigma_m (black) and sqrt(sum sigma_c^2) open circle — the
        gap visualises the screening cancellation.
Right: decomposed PSD (mean over 8 sites) for the 4 sources + measured total
        (grey). Annotate integral ratio sum(PSD_comp)/PSD_total ~= 2.0.
"""
from __future__ import annotations
import sys, csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (setup, save, TRP_NAMES, SOURCE_COLORS, SOURCE_LABELS,
                     PAPER_DIR)
from utils import trapz

P1 = PAPER_DIR.parent / 'phase1_basic_stats'
P2 = PAPER_DIR.parent / 'phase2_timescale'
P3 = PAPER_DIR.parent / 'phase3_source_attribution'


def load_sigma():
    sig = np.zeros((8, 5))   # cols: total, protein, water, nucleotide, ions
    with open(P1 / 'sigma_matrix.csv') as f:
        for r in csv.DictReader(f):
            if r['dataset'] != 'slow':
                continue
            i = int(r['idx'])
            sig[i] = [float(r['sigma_total_cm']), float(r['sigma_protein_cm']),
                      float(r['sigma_water_cm']), float(r['sigma_nucleotide_cm']),
                      float(r['sigma_ions_cm'])]
    return sig


def main():
    setup()
    sig = load_sigma()
    src_rms = sig[:, 1:]                       # (8,4)
    sigma_tot = sig[:, 0]
    sigma_unc = np.sqrt((src_rms**2).sum(axis=1))

    # PSDs
    comp = np.load(P3 / 'comp_psd_stitched.npz')
    f = comp['f_cm']
    comp_psd = {s: comp[s].mean(axis=1) for s in SOURCE_LABELS}   # mean over sites
    tot = np.load(P2 / 'psd_stitched.npz')
    psd_tot = tot['psd'].mean(axis=1)
    sum_comp = sum(comp_psd.values())
    ratio = trapz(sum_comp, f) / trapz(psd_tot, f)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6),
                                   gridspec_kw={'width_ratios': [1.05, 1]})

    # ── left: grouped bars + screening markers ──
    x = np.arange(8); w = 0.18
    for ig, s in enumerate(SOURCE_LABELS):
        ax1.bar(x + (ig - 1.5) * w, src_rms[:, ig], w,
                color=SOURCE_COLORS[s], edgecolor='black', linewidth=0.3,
                label=s.capitalize())
    ax1.scatter(x, sigma_tot, color='black', marker='_', s=150, zorder=5,
                linewidth=2.2, label=r'total $\sigma_m$')
    ax1.scatter(x, sigma_unc, facecolor='none', edgecolor='#444', marker='o',
                s=75, zorder=4, linewidth=1.3,
                label=r'$\sqrt{\sum_g\sigma_{g,m}^{2}}$ (uncorrelated)')
    # join each total to its uncorrelated circle to emphasise the gap
    for i in range(8):
        ax1.plot([x[i], x[i]], [sigma_tot[i], sigma_unc[i]],
                 color='#444', lw=0.7, alpha=0.5, zorder=3)
    ax1.set_xticks(x); ax1.set_xticklabels(TRP_NAMES)
    ax1.set_ylabel(r'RMS $\delta\epsilon_m$ (cm$^{-1}$)')
    ax1.set_title('(a) Per-source contribution & screening')
    ax1.set_ylim(0, 3200)
    ax1.legend(loc='upper left', ncol=3, columnspacing=1, handlelength=1.2)
    ax1.grid(axis='y', alpha=0.25)
    # RMS-based uncorrelated-sum / total ratio (mirrors the right-panel PSD ratio)
    rms_ratio = float(np.mean(sigma_unc / sigma_tot))
    ax1.text(0.97, 0.90,
             r'$\langle\sqrt{\sum_g\sigma_{g,m}^{2}}\,/\,\sigma_m\rangle = $'
             f'{rms_ratio:.2f}',
             transform=ax1.transAxes, va='top', ha='right', fontsize=9,
             bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#bbb', alpha=0.9))

    # ── right: decomposed PSD ──
    # black = measured total; grey = uncorrelated sum of components (screening gap)
    ax2.loglog(f, psd_tot, '-', color='black', lw=2.0, alpha=0.95,
               label='total (measured)')
    ax2.loglog(f, sum_comp, '--', color='grey', lw=1.5, alpha=0.9,
               label=r'$\sum_g S_g$ (uncorrelated)')
    for s in SOURCE_LABELS:
        ax2.loglog(f, comp_psd[s], '-', color=SOURCE_COLORS[s], lw=1.1,
                   alpha=0.8, label=s.capitalize())
    ax2.axvline(1.0, color='gray', ls=':', lw=0.8, alpha=0.6)
    ax2.set_xlabel(r'frequency $f$ (cm$^{-1}$)')
    ax2.set_ylabel('PSD (a.u.)')
    ax2.set_title('(b) Source-resolved PSD (mean over sites)')
    ax2.legend(loc='lower left', ncol=2, columnspacing=1, handlelength=1.2)
    ax2.text(0.97, 0.97,
             r'$\int\sum_g S_g\,df\;/\;\int S_{\rm tot}\,df = $'
             f'{ratio:.2f}',
             transform=ax2.transAxes, va='top', ha='right', fontsize=9,
             bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#bbb', alpha=0.9))

    fig.tight_layout()
    save(fig, 'fig2_source_decomp.png')
    # screening numbers for the writeup
    Rscreen = 1 - sigma_tot**2 / (src_rms**2).sum(axis=1)
    print(f"R_screen mean = {Rscreen.mean():.3f}  (range {Rscreen.min():.2f}-{Rscreen.max():.2f})")
    print(f"sigma_uncorr/sigma_tot mean = {(sigma_unc/sigma_tot).mean():.3f}")
    print(f"PSD integral ratio = {ratio:.3f}")


if __name__ == '__main__':
    main()
