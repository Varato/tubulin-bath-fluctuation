#!/usr/bin/env python
"""Fig 1: stitched ACF (with tri-exp fit) + stitched PSD, all 8 Trps.

Left:  ACF, dots = raw stitched, solid = tri-exp model from CSV params. x = ps (log).
Right: PSD, lines. x = cm^-1 (log), y = a.u. (log).
"""
from __future__ import annotations
import sys, csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (setup, save, TRP_NAMES, TRP_COLORS, PAPER_DIR)

P2 = PAPER_DIR.parent / 'phase2_timescale'
T_CROSS_PS = 10.0   # stitch crossover


def load_triexp():
    """Corrected tri-exp (tau3 anchored to slow-only) — same source as Table 1."""
    amps = np.zeros((8, 3)); taus = np.zeros((8, 3))
    with open(PAPER_DIR / 'corrected_triexp.csv') as f:
        for i, r in enumerate(csv.DictReader(f)):
            amps[i] = [float(r['A1']), float(r['A2']), float(r['A3'])]
            taus[i] = [float(r['tau1_ps']), float(r['tau2_ps']), float(r['tau3_ps'])]
    return amps, taus


def main():
    setup()
    acf_d = np.load(P2 / 'acf_stitched.npz')
    t_ps = acf_d['t_lag_ps']; acf = acf_d['acf']          # (N,8)
    psd_d = np.load(P2 / 'psd_stitched.npz')
    f_cm = psd_d['f_cm']; psd = psd_d['psd']              # (N,8)
    amps, taus = load_triexp()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

    # ── left: ACF + fit ──
    for i in range(8):
        c = TRP_COLORS[i]
        ax1.plot(t_ps, acf[:, i], 'o', color=c, ms=2.8, alpha=0.45,
                 markeredgewidth=0)
        fit = sum(amps[i, k] * np.exp(-t_ps / taus[i, k]) for k in range(3))
        ax1.plot(t_ps, fit, '-', color=c, lw=1.4, alpha=0.95, label=TRP_NAMES[i])
    ax1.axvline(T_CROSS_PS, color='gray', ls=':', lw=0.8, alpha=0.6)
    ax1.text(T_CROSS_PS*1.15, 0.5, '10 ps\nstitch', fontsize=7, color='gray')
    ax1.set_xscale('log'); ax1.set_ylim(-0.05, 1.05)
    ax1.set_xlabel('lag time $t$ (ps)')
    ax1.set_ylabel(r'$C(t)\,/\,C(0)$')
    ax1.set_title('(a) Autocorrelation')
    ax1.legend(ncol=2, loc='upper right', columnspacing=1, handlelength=1.2)

    # ── right: PSD ──
    for i in range(8):
        ax2.plot(f_cm, psd[:, i], '-', color=TRP_COLORS[i], lw=1.2,
                 alpha=0.85, label=TRP_NAMES[i])
    ax2.axvline(1.0, color='gray', ls=':', lw=0.8, alpha=0.6)
    ax2.text(1.1, ax2.get_ylim()[1]*0.3 if ax2.get_yscale()=='linear' else 1e7,
             '1 cm⁻¹\nstitch', fontsize=7, color='gray')
    ax2.set_xscale('log'); ax2.set_yscale('log')
    ax2.set_xlabel(r'frequency $f$ (cm$^{-1}$)')
    ax2.set_ylabel('PSD (a.u.)')
    ax2.set_title('(b) Power spectral density')

    fig.tight_layout()
    save(fig, 'fig1_acf_psd.png')


if __name__ == '__main__':
    main()
