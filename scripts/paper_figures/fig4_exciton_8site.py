#!/usr/bin/env python
"""Fig 4: 8-site exciton dynamics on the full Craddock Trp network.

Left:  coloured-noise trajectory MC (full tri-exp OU bath per site). All 8 P_m(t)
       shown; Trp4 (donor) and Trp7 (strongest-coupled neighbour, J=-59 cm^-1)
       emphasised in black/red with SEM bands.
Right: Haken-Strobl Lindblad (gamma=50 cm^-1, Craddock 2014) on the same
       H_CRAD. Same colour scheme, dashed.

Direct visual comparison: under MD noise, excitation stays concentrated on
Trp4 / Trp7 (the dimer subsystem); under Markovian Lindblad, the bath
over-delocalises the population across the network.

Annotates P_leak(2 ps) = 1 - P_Trp4 - P_Trp7 in each panel: the fraction of
excitation that has escaped the Trp4-Trp7 dimer into the other six sites.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (setup, save, TRP_NAMES, TRP_COLORS, PAPER_DIR)


def main():
    setup()
    mc = np.load(PAPER_DIR / 'exciton_mc_8site.npz')
    lb = np.load(PAPER_DIR / 'lindblad_8site.npz')
    t = mc['t']; t_fs = t * 1e3
    P_mc = mc['P_full']           # (8, N_STEPS)
    P_mc_sem = mc['P_full_sem']
    P_lb = lb['P']                # (8, N_STEPS)
    donor = int(mc['donor'])      # 3 = Trp4
    target = int(mc['target'])    # 6 = Trp7

    # ── muted treatment for the 6 non-focal sites ──
    FOCAL = {donor, target}
    focal_style = {donor: ('black', 'Trp4 (donor)'),
                   target: ('#d62728', 'Trp7 (target)')}

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(11, 4.6), sharey=True,
        gridspec_kw={'width_ratios': [1, 1], 'wspace': 0.08})
    ax2.tick_params(labelleft=False)

    def plot_panel(ax, P, sem=None, title='', dashed=False):
        # non-focal sites first (muted, thin, no band)
        for m in range(8):
            if m in FOCAL:
                continue
            ls = '--' if dashed else '-'
            ax.plot(t_fs, P[m], ls, color=TRP_COLORS[m], lw=1.0, alpha=0.75,
                    label=TRP_NAMES[m])
        # focal sites on top (heavy, with SEM band on MC)
        for m, (col, lbl) in focal_style.items():
            ls = '--' if dashed else '-'
            ax.plot(t_fs, P[m], ls, color=col, lw=2.0, alpha=0.95, label=lbl)
            if sem is not None and not dashed:
                ax.fill_between(t_fs, P[m] - sem[m], P[m] + sem[m],
                                color=col, alpha=0.15)
        ax.set_xlabel('time (fs)')
        ax.set_xlim(0, t_fs[-1])
        ax.set_ylim(-0.02, 1.0)
        ax.set_title(title)
        # P_leak annotation
        P_dimer_end = P[donor, -1] + P[target, -1]
        ax.text(0.03, 0.97,
                f'$P_{{\\rm leak}}(2\\,$ps$) = {1 - P_dimer_end:.2f}$',
                transform=ax.transAxes, va='top', ha='left', fontsize=12,
                bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#bbb',
                          alpha=0.9))

    plot_panel(ax1, P_mc, sem=P_mc_sem,
               title='(a) Coloured-noise MC (MD bath)')
    ax1.set_ylabel(r'population $P_m(t)$')
    ax1.legend(loc='center right', fontsize=10)

    plot_panel(ax2, P_lb, dashed=True,
               title=r'(b) Haken-Strobl Lindblad, $\gamma=50$ cm$^{-1}$')
    ax2.legend(loc='center right', fontsize=10)

    # ── final-distribution summary printed for the writeup ──
    print("P_m(2 ps):")
    print(f"  {'site':<6} {'MC':>8} {'Lindblad':>10} {'diff':>8}")
    for m in range(8):
        print(f"  {TRP_NAMES[m]:<6} {P_mc[m,-1]:>8.3f} {P_lb[m,-1]:>10.3f} "
              f"{P_mc[m,-1]-P_lb[m,-1]:>+8.3f}")
    print(f"  P_leak(2ps):  MC={1-P_mc[donor,-1]-P_mc[target,-1]:.3f}  "
          f"Lindblad={1-P_lb[donor,-1]-P_lb[target,-1]:.3f}")

    fig.tight_layout()
    save(fig, 'fig4_exciton_8site.png')


if __name__ == '__main__':
    main()
