#!/usr/bin/env python
"""Fig A2: site-energy fluctuation traces (4x2, one panel per Trp).

Each panel overlays the last 2 ns of the slow trajectory (coarse, 10 ps steps)
with the full 2 ns fast trajectory (10 fs steps, downsampled) on a common
2 ns time axis. Both are shown in cm^-1, demeaned.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import setup, save, TRP_NAMES, TRP_COLORS, PAPER_DIR
from utils import load, V_TO_CM


def main():
    setup()
    slow = load('slow')
    fast = load('fast')
    ds_slow = slow['delta_s_total'] * V_TO_CM          # (Ns, 8)
    ds_fast = fast['delta_s_total'] * V_TO_CM          # (Nf, 8)
    dt_slow = slow['dt_ps']; dt_fast = fast['dt_ps']

    # last 2 ns of slow = last 2000 ps / 10 ps = last 200 frames
    n_win = int(2000.0 / dt_slow)
    sl = ds_slow[-n_win:]
    t_slow = np.arange(n_win) * dt_slow - n_win * dt_slow   # shift to end at 0... keep as 0..2ns
    t_slow = np.arange(n_win) * dt_slow                      # 0..2000 ps

    # fast full 2 ns, downsample 10 fs -> 1 ps (every 100th) for plotting
    step = max(1, int(1.0 / dt_fast))
    ff = ds_fast[::step]
    t_fast = np.arange(len(ff)) * dt_fast * step

    fig, axs = plt.subplots(4, 2, figsize=(11, 9.5), sharex=True)
    for i in range(8):
        ax = axs[i // 2, i % 2]
        ax.plot(t_fast, ff[:, i] - ff[:, i].mean(), color=TRP_COLORS[i],
                lw=0.6, alpha=0.55, label='fast (10 fs)')
        ax.plot(t_slow, sl[:, i] - sl[:, i].mean(), color='black',
                lw=1.3, alpha=0.9, label='slow (10 ps)')
        ax.set_ylabel(TRP_NAMES[i], fontsize=13, rotation=0, labelpad=22,
                      va='center')
        if i == 0:
            ax.legend(loc='upper right')
        if i >= 6:
            ax.set_xlabel('window time (ps)')
        ax.set_xlim(0, 2000)
    fig.suptitle('Site-energy fluctuation traces (2 ns window, demeaned)',
                 y=0.995, fontsize=16)
    fig.tight_layout()
    save(fig, 'figA2_traces.png')
    print(f"slow window: {n_win} frames; fast plotted: {len(ff)} frames")


if __name__ == '__main__':
    main()
