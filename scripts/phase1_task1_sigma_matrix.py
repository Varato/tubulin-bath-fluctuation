#!/usr/bin/env python
"""Phase 1 / Task 1.1 — fluctuation amplitude & disorder strength (RESEARCH_PLAN §2.2).

Per-site σ for total + 4 component-resolved site-energy fluctuations on both
datasets (cm^-1). Disorder-coupling ratio σ_total / J (J = 45 cm^-1) confirms
the system is in the σ/J ≫ 1 strong-disorder regime.

Outputs (results/phase1_basic_stats/):
  sigma_matrix.csv     dataset,site,idx,{total,protein,water,nucleotide,ions} in cm^-1
  disorder_ratio.csv   dataset,site,idx,sigma_total_cm,sigma_over_J
  sigma_matrix.png     σ_total + component breakdown per site; σ/J panel
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load, V_TO_CM, J_COUPLING, TRP_PRETTY, COMPONENTS,
                   phase_dir, write_csv, setup_style, save_fig)

COMPS = ['total'] + COMPONENTS   # 5 columns


def sigma_table(tag):
    d = load(tag)
    # shape (N, 8) per component; std over time axis 0 → (8,)
    return {c: d[f'delta_s_{c}'].std(axis=0) * V_TO_CM for c in COMPS}, d['n_frames']


rows_mat, rows_ratio = [], []
sig = {t: {} for t in ('slow', 'fast')}
for tag in ('slow', 'fast'):
    tbl, N = sigma_table(tag)
    sig[tag] = tbl
    for s in range(8):
        rec = [tag, TRP_PRETTY[s], s] + [f"{tbl[c][s]:.3f}" for c in COMPS]
        rows_mat.append(rec)
        rows_ratio.append((tag, TRP_PRETTY[s], s,
                           f"{tbl['total'][s]:.3f}",
                           f"{tbl['total'][s] / J_COUPLING:.3f}"))
    # system summary
    sig_tot = tbl['total']
    print(f"[{tag}] N={N}")
    print(f"   σ_total (cm^-1): per-site min={sig_tot.min():.1f} "
          f"mean={sig_tot.mean():.1f}±{sig_tot.std(ddof=0):.1f} max={sig_tot.max():.1f}")
    print(f"   σ_total / J  :   min={sig_tot.min()/J_COUPLING:.2f} "
          f"mean={sig_tot.mean()/J_COUPLING:.2f}±{sig_tot.std(ddof=0)/J_COUPLING:.2f} "
          f"max={sig_tot.max()/J_COUPLING:.2f}   (J={J_COUPLING} cm^-1)")
    # component contributions to σ (not σ² — just for orientation)
    print("   component σ means (cm^-1): " +
          "  ".join(f"{c}={tbl[c].mean():.1f}" for c in COMPONENTS))
    print()

P = phase_dir(1)
write_csv(P / 'sigma_matrix.csv',
          "dataset,site,idx," + ",".join(f"sigma_{c}_cm" for c in COMPS),
          rows_mat)
write_csv(P / 'disorder_ratio.csv',
          "dataset,site,idx,sigma_total_cm,sigma_over_J", rows_ratio)
print(f"wrote {P/'sigma_matrix.csv'} and {P/'disorder_ratio.csv'}")

# ------------------------------------------------------------------ figure
setup_style()
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 6.2), sharex=True)
xs = np.arange(8)
w = 0.38
# panel 1: σ_total slow vs fast
ax1.bar(xs - w/2, sig['slow']['total'], w, label=f'slow (50 ns, N={load("slow")["n_frames"]})', color='#225ea8')
ax1.bar(xs + w/2, sig['fast']['total'], w, label=f'fast (2 ns,  N={load("fast")["n_frames"]})', color='#41b6c4')
ax1.axhline(J_COUPLING, color='crimson', ls='--', lw=1, label=f'J = {J_COUPLING} cm⁻¹')
ax1.set_ylabel(r'$\sigma_\mathrm{total}$ (cm⁻¹)')
ax1.set_title('Phase 1 / Task 1.1 — disorder strength per Trp site')
ax1.legend(loc='upper right', fontsize=9)

# panel 2: σ/J on log scale to show strong-disorder regime
for off, tag, col in ((-w/2, 'slow', '#225ea8'), (w/2, 'fast', '#41b6c4')):
    r = sig[tag]['total'] / J_COUPLING
    ax2.bar(xs + off, r, w, color=col)
ax2.axhline(1.0, color='crimson', ls='--', lw=1, label='σ/J = 1 (weak→strong crossover)')
ax2.set_yscale('log')
ax2.set_ylabel(r'$\sigma_\mathrm{total}\,/\,J$')
ax2.set_xticks(xs)
ax2.set_xticklabels(TRP_PRETTY, rotation=20, ha='right')
ax2.legend(loc='upper right', fontsize=9)

save_fig(fig, P / 'sigma_matrix.png')
print(f"wrote {P/'sigma_matrix.png'}")
