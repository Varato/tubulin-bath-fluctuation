#!/usr/bin/env python
"""Phase 3 / Task 3.1 supplement — tri-exp per component figure.

Regenerates the component ACF figure showing the TRI-EXP decomposition per
component (the analysis that maps to Phase 2's τ₁/τ₂/τ₃), instead of the
single-exp τ_fast/τ_slow shown in comp_acf.png.

Each panel: one component's stitched ACF (mean over sites) + the three
tri-exp components colour-coded by the Phase 2 correspondence.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (COMPONENTS, TRP_PRETTY, phase_dir, setup_style, save_fig,
                   log_interp)

T_CROSS = 10.0


def triexp(t, A1, A2, t1, t2, t3):
    return A1*np.exp(-t/t1) + A2*np.exp(-t/t2) + (1-A1-A2)*np.exp(-t/t3)


def main():
    P = phase_dir(3)
    fast = np.load(P / 'comp_acf_fast.npz')
    slow = np.load(P / 'comp_acf_slow.npz')

    t_log = np.logspace(np.log10(0.01), np.log10(10000), 250)
    m = t_log > 0

    # fit + collect curves
    fits = {}
    for comp in COMPONENTS:
        tf = fast[f't_{comp}']; af = fast[comp]
        ts = slow[f't_{comp}']; asl = slow[comp]
        a_fast = log_interp(t_log, tf, af.mean(axis=1))
        a_slow = log_interp(t_log, ts, asl.mean(axis=1))
        stitch = np.where(t_log < T_CROSS, a_fast, a_slow)
        try:
            p, _ = curve_fit(triexp, t_log[m], stitch[m],
                             p0=[0.4, 0.3, 0.05, 1.5, 800],
                             bounds=([0, 0, 1e-3, 0.1, 5], [1, 1, 5, 100, 5e4]),
                             maxfev=80000)
            pairs = sorted([(p[0], p[2]), (p[1], p[3]),
                            (1 - p[0] - p[1], p[4])], key=lambda x: x[1])
            fits[comp] = dict(stitch=stitch, pairs=pairs, p=p)
        except Exception:
            fits[comp] = dict(stitch=stitch, pairs=None, p=None)

    setup_style()
    import matplotlib.pyplot as plt

    colors_comp = {'protein': '#6a51a3', 'water': '#41b6c4',
                   'nucleotide': '#fd8d3c', 'ions': '#74c476'}
    # component colours for the three modes, keyed by Phase 2 correspondence
    c_lib = '#fd8d3c'    # τ₁ libration (orange)
    c_rot = '#74c476'    # τ₂ rotation  (green)
    c_slow = '#6a51a3'   # τ₃ slow       (purple)

    fig, axs = plt.subplots(2, 2, figsize=(10, 6.8), sharex=True)
    for k, comp in enumerate(COMPONENTS):
        ax = axs[k // 2, k % 2]
        st = fits[comp]['stitch']
        ax.semilogx(t_log, st, color=colors_comp[comp], lw=1.4, alpha=0.6,
                    label='stitched ACF')
        pairs = fits[comp]['pairs']
        if pairs is not None:
            total = triexp(t_log, *fits[comp]['p'])
            ax.semilogx(t_log, total, 'r-', lw=1.0, alpha=0.8, label='tri-exp fit')
            labels = [r'$\tau_1$ (libration)', r'$\tau_2$ (rotation)',
                      r'$\tau_3$ (slow)']
            cc = [c_lib, c_rot, c_slow]
            for i, (A, tau) in enumerate(pairs):
                ax.semilogx(t_log, A * np.exp(-t_log / tau), color=cc[i],
                            lw=0.9, ls='--',
                            label=fr'{labels[i]}: $\tau$={tau:.3f} ps, A={A:.2f}')
            title = f"{comp}   τ₁={pairs[0][1]:.3f}  τ₂={pairs[1][1]:.2f}  τ₃={pairs[2][1]:.0f} ps"
        else:
            title = f"{comp}  (fit failed)"
        ax.axvline(T_CROSS, color='gray', ls=':', lw=0.6)
        ax.set_title(title, fontsize=8.5)
        ax.set_ylim(-0.1, 1.05)
        if k == 0:
            ax.legend(fontsize=6, loc='upper right')
        if k // 2 == 1:
            ax.set_xlabel('lag t (ps)')
        if k % 2 == 0:
            ax.set_ylabel('C(t)')

    fig.suptitle(
        'Phase 3 / Task 3.1 — tri-exp per component (stitched ACF, mean over 8 sites)\n'
        'dashed = three components;  grey : = 10 ps crossover\n'
        'Phase 2 reference: τ₁=0.044 (libration)  τ₂=1.70 (rotation)  τ₃=1140–2663 ps (protein)',
        y=1.0, fontsize=9)
    save_fig(fig, P / 'comp_acf_triexp.png')
    print(f"wrote {P/'comp_acf_triexp.png'}")

    # also: a clean correspondence summary bar chart
    fig2, ax = plt.subplots(figsize=(7, 4))
    p2 = [0.044, 1.70, 2663]
    labels_p2 = ['τ₁ libration\n(0.044)', 'τ₂ rotation\n(1.70)', 'τ₃ protein\n(2663)']
    water_taus = [fits['water']['pairs'][i][1] if fits['water']['pairs'] else np.nan
                  for i in range(3)]
    prot_taus = [fits['protein']['pairs'][i][1] if fits['protein']['pairs'] else np.nan
                 for i in range(3)]
    x = np.arange(3)
    w = 0.25
    ax.bar(x - w, p2, w, color='k', alpha=0.5, label='Phase 2 (total ACF)')
    ax.bar(x, water_taus, w, color='#41b6c4', label='water (per-component)')
    ax.bar(x + w, prot_taus, w, color='#6a51a3', label='protein (per-component)')
    ax.set_yscale('log')
    ax.set_xticks(x); ax.set_xticklabels(labels_p2)
    ax.set_ylabel('τ (ps, log scale)')
    ax.set_title('Phase 2 ↔ Phase 3 timescale correspondence')
    ax.legend(fontsize=8)
    save_fig(fig2, P / 'tau_correspondence.png')
    print(f"wrote {P/'tau_correspondence.png'}")


if __name__ == '__main__':
    main()
