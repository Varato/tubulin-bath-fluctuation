#!/usr/bin/env python
"""Fig 3: exciton population dynamics in the Trp4-Trp7 dimer.

(a) Component ablation: P7(t) for noiseless, tau1, full, tau3-only. SEM bands.
    Shows fast noise enables incoherent transfer (~0.32) while static disorder
    (tau3) traps the excitation (~0.06).
(b) Noise-model comparison: trajectory MC (full, colored noise from MD) vs the
    Craddock 2014 Haken-Strobl Lindblad (gamma=50 cm^-1) on the same dimer.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import setup, save, PAPER_DIR


def main():
    setup()
    mc = np.load(PAPER_DIR / 'exciton_mc.npz')
    lb2 = np.load(PAPER_DIR / 'lindblad_2site.npz')
    t = mc['t']; t_fs = t * 1e3

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.6), sharey=True,
                                   gridspec_kw={'width_ratios': [1, 1],
                                                'wspace': 0.08})
    ax2.tick_params(labelleft=False)   # shared y-axis: drop duplicate labels

    # ── (a) ablation ──
    ax1.plot(t_fs, mc['P7_noiseless'], ':', color='gray', lw=1.6, alpha=0.7,
             label='noiseless')
    for tag, col, lbl in [
        ('tau1', '#d62728', r'$\tau_1$ only (libration)'),
        ('full', 'black',   'full ($\\tau_1{+}\\tau_2{+}\\tau_3$)'),
        ('tau3', '#2ca02c', r'$\tau_3$ only (static)')]:
        P7 = mc[f'{tag}_P7']; sem = mc[f'{tag}_P7_sem']
        ax1.plot(t_fs, P7, '-', color=col, lw=1.5, label=lbl)
        ax1.fill_between(t_fs, P7 - sem, P7 + sem, color=col, alpha=0.15)
    ax1.set_xlabel('time (fs)'); ax1.set_ylabel(r'population $P_7(t)$')
    ax1.set_title('(a) Noise-component ablation')
    ax1.set_ylim(-0.02, 1.0); ax1.set_xlim(0, t_fs[-1])
    ax1.legend(loc='center right')
    # annotate plateau values at 2 ps
    txt = (f"$P_7(2\\,$ps$)$:  full={mc['full_P7'][-1]:.2f},  "
           f"$\\tau_3$={mc['tau3_P7'][-1]:.2f}")
    ax1.text(0.03, 0.97, txt, transform=ax1.transAxes, va='top', fontsize=8.5,
             bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#bbb', alpha=0.9))

    # ── (b) MC vs Lindblad ──
    ax2.plot(t_fs, mc['full_P7'], '-', color='black', lw=1.7,
             label='trajectory MC (MD noise)')
    ax2.fill_between(t_fs, mc['full_P7'] - mc['full_P7_sem'],
                     mc['full_P7'] + mc['full_P7_sem'], color='black', alpha=0.15)
    ax2.plot(lb2['t'] * 1e3, lb2['P7'], '--', color='#1f77b4', lw=1.7,
             label=r'Lindblad, $\gamma=50$ cm$^{-1}$')
    ax2.axhline(0.5, color='#1f77b4', ls=':', lw=0.7, alpha=0.5)
    ax2.set_xlabel('time (fs)')
    ax2.set_title('(b) Exact MD noise vs Markovian Lindblad')
    ax2.set_xlim(0, t_fs[-1])
    ax2.legend(loc='lower right')
    txt2 = (f"$P_7(2\\,$ps$)$:  MC={mc['full_P7'][-1]:.2f},  "
            f"Lindblad={lb2['P7'][-1]:.2f}")
    ax2.text(0.03, 0.97, txt2, transform=ax2.transAxes, va='top', fontsize=8.5,
             bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#bbb', alpha=0.9))

    fig.tight_layout()
    save(fig, 'fig3_exciton.png')
    print(f"noiseless peak = {mc['P7_noiseless'].max():.3f}")
    print(f"P7(2ps): tau1={mc['tau1_P7'][-1]:.3f} tau2={mc['tau2_P7'][-1]:.3f} "
          f"tau1+2={mc['tau1_tau2_P7'][-1]:.3f} full={mc['full_P7'][-1]:.3f} "
          f"tau3={mc['tau3_P7'][-1]:.3f}")
    print(f"Lindblad 2-site P7(2ps)={lb2['P7'][-1]:.3f}")


if __name__ == '__main__':
    main()
