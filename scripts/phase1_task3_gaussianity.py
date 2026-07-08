#!/usr/bin/env python
"""Phase 1 / Task 1.3 — Gaussianity test (RESEARCH_PLAN §2.4).

If site-energy fluctuations are Gaussian, the ACF/PSD used in Phase 2 are
complete second-order descriptors. Tests this on the fast dataset (full
sub-ps resolution).

Outputs (results/phase1_basic_stats/):
  gaussianity.csv    per-site skewness, excess kurtosis, KS D & p-value
  gaussianity.png    2×4 histograms with Gaussian overlay
  gaussian_zscore.npz  z-scored series (handy for Phase 5 sampling)
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load, V_TO_CM, TRP_PRETTY, phase_dir, write_csv,
                   setup_style, save_fig)

d = load('fast')
x = d['delta_s_total'] * V_TO_CM    # (N, 8) cm^-1
N = d['n_frames']

rows, z = [], np.empty_like(x)
for s in range(8):
    xs = x[:, s]
    mu, sigma = xs.mean(), xs.std()
    z[:, s] = (xs - mu) / sigma
    skew = stats.skew(xs)
    kurt = stats.kurtosis(xs, fisher=True)   # excess kurtosis (Gaussian → 0)
    # KS vs fitted normal (note: params estimated → Lilliefors regime, p biased low)
    D, p = stats.kstest(z[:, s], 'norm')
    rows.append((TRP_PRETTY[s], s, f"{skew:+.4f}", f"{kurt:+.4f}",
                 f"{D:.4f}", f"{p:.2e}"))
    print(f"  {TRP_PRETTY[s]:7s}  skew={skew:+.4f}  "
          f"excess_kurt={kurt:+.4f}  KS D={D:.4f}  p={p:.2e}")

print(f"\nN = {N} frames. Caveat: with N={N:.1e}, formal tests reject Gaussianity "
      f"for any tiny deviation — judge by the magnitudes (skew, kurt).")

P = phase_dir(1)
write_csv(P / 'gaussianity.csv',
          "site,idx,skewness,excess_kurtosis,ks_D,ks_pvalue", rows)
np.savez_compressed(P / 'gaussian_zscore.npz',
                    z=z, time_ps=d['time_ps'], dt_ps=d['dt_ps'],
                    labels=np.array(TRP_PRETTY))
print(f"wrote {P/'gaussianity.csv'} and {P/'gaussian_zscore.npz'}")

# ------------------------------------------------------------------ figure
setup_style()
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 4, figsize=(11, 5.2), sharey=False)
xs_grid = np.linspace(-4, 4, 400)
for s in range(8):
    ax = axs[s // 4, s % 4]
    ax.hist(z[:, s], bins=120, density=True, color='#41b6c4', alpha=0.65)
    ax.plot(xs_grid, stats.norm.pdf(xs_grid), color='crimson', lw=1.4, label='N(0,1)')
    ax.set_title(f"{TRP_PRETTY[s]}  | skew {rows[s][2]}  kurt {rows[s][3]}",
                 fontsize=9)
    ax.set_xlim(-4, 4)
    ax.legend(fontsize=8, loc='upper right')
    if s % 4 == 0:
        ax.set_ylabel('density')
    if s // 4 == 1:
        ax.set_xlabel('z-score')
fig.suptitle('Phase 1 / Task 1.3 — Gaussianity check (fast dataset, '
             f'N={N})', y=1.00)
save_fig(fig, P / 'gaussianity.png')
print(f"wrote {P/'gaussianity.png'}")
