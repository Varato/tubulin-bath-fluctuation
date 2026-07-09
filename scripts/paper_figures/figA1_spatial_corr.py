#!/usr/bin/env python
"""Fig A1: 8x8 spatial correlation matrix of site-energy fluctuations (slow traj).

Mean off-diagonal |r| ~ 0.02 -> sites are independent.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import setup, save, TRP_NAMES, PAPER_DIR

P4 = PAPER_DIR.parent / 'phase4_spatial_correlation'


def main():
    setup()
    d = np.load(P4 / 'corr_cov_slow.npz')
    # 'total' stored as covariance; derive correlation. Check if a pearson exists.
    key = 'total'
    cov = d[key]
    sig = np.sqrt(np.diag(cov))
    corr = cov / np.outer(sig, sig)

    off = corr.copy()
    np.fill_diagonal(off, np.nan)
    mean_off = np.nanmean(np.abs(off))

    fig, ax = plt.subplots(figsize=(6.6, 6.0))
    im = ax.imshow(corr, cmap='RdBu_r', vmin=-0.4, vmax=0.4)
    ax.set_xticks(range(8)); ax.set_xticklabels(TRP_NAMES)
    ax.set_yticks(range(8)); ax.set_yticklabels(TRP_NAMES)
    for i in range(8):
        for j in range(8):
            ax.text(j, i, f'{corr[i,j]:.2f}', ha='center', va='center',
                    fontsize=12, color='black')
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label('Pearson $r$')
    ax.set_title(f'Spatial correlation (slow traj)\n'
                 f'mean off-diag $|r|$ = {mean_off:.3f}')
    fig.tight_layout()
    save(fig, 'figA1_spatial_corr.png')
    print(f"mean off-diag |r| = {mean_off:.4f}")


if __name__ == '__main__':
    main()
