#!/usr/bin/env python
"""4×2 figure: MC vs Lindblad from each starting Trp on the 8-site network.

Loads precomputed data from run_8site_scan.py.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import setup, PAPER_DIR

N_SITE = 8
COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
LABELS = [f'Trp{i+1}' for i in range(N_SITE)]


def main():
    setup()

    d = np.load(PAPER_DIR / 'fig4_8site_scan.npz')
    t_fs = d['t'] * 1e3
    mc = d['mc']   # (8, 8, N_STEPS)
    lb = d['lb']

    fig, axes = plt.subplots(4, 2, figsize=(10, 10), sharex=True, sharey=True)
    plt.subplots_adjust(left=0.08, right=0.96, top=0.94, bottom=0.06,
                        hspace=0.18, wspace=0.10)

    for s in range(N_SITE):
        ax = axes[s // 2, s % 2]
        for n in range(N_SITE):
            ax.plot(t_fs, mc[s, n], '-', color=COLORS[n], lw=1.8,
                    alpha=0.9 if n == s else 0.55)
            ax.plot(t_fs, lb[s, n], '--', color=COLORS[n], lw=1.6,
                    alpha=0.8 if n == s else 0.45)
        ax.set_title(f'Initial: Trp{s+1}', fontsize=12)
        if s // 2 == 3:
            ax.set_xlabel('time (fs)')
        if s % 2 == 0:
            ax.set_ylabel('population')
        ax.set_xlim(0, t_fs[-1])
        ax.set_ylim(-0.02, 1.0)

    legend_lines = ([Line2D([0], [0], color=COLORS[n], lw=2.0,
                            label=LABELS[n]) for n in range(N_SITE)]
                    + [Line2D([0], [0], color='gray', lw=2.0, ls='-',
                              label='MC (colored noise)'),
                       Line2D([0], [0], color='gray', lw=2.0, ls='--',
                              label='Lindblad ($\\gamma=50$)')])
    axes[0, 0].legend(handles=legend_lines, ncol=3, fontsize=7.5,
                      loc='center left')

    p_png = PAPER_DIR / 'fig4_8site_scan.png'
    p_pdf = PAPER_DIR / 'fig4_8site_scan.pdf'
    fig.savefig(p_png, dpi=300)
    fig.savefig(p_pdf)
    plt.close(fig)
    print(f"saved {p_png}")


if __name__ == '__main__':
    main()
