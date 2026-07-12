#!/usr/bin/env python
"""2-panel figure: (a) clean eigenstate PRs, (b) disorder sigma-scan of max PR.

Replaces fig_oscillator_strength.py. Uses only eigenvector data (no dipole
directions needed), so the analysis is exact.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.use('Agg')
plt.rcParams.update({
    'figure.dpi': 130, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'font.size': 16, 'axes.labelsize': 18, 'axes.titlesize': 17,
    'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 12,
    'axes.grid': True, 'grid.alpha': 0.25, 'axes.axisbelow': True,
    'lines.linewidth': 1.5, 'legend.frameon': False,
    'font.family': 'DejaVu Sans',
})

RESULTS_DIR = Path(__file__).resolve().parent.parent / 'results'
FIG_DIR = Path(__file__).resolve().parent.parent / 'figures'

# Load data
clean = np.load(RESULTS_DIR / 'clean_oscillator_strength.npz')
scan = np.load(RESULTS_DIR / 'disorder_scan.npz')

energies = clean['energies']
pr_clean = clean['pr']

sigmas = scan['sigmas']
pr_mean = scan['pr_mean']  # shape (n_sigmas, 8)
pr_sem = scan['pr_sem']

SIGMA_3 = 200.0  # cm^-1

# Max PR across eigenstates at each sigma
pr_max_mean = np.max(pr_mean, axis=1)
# SEM: propagate from the eigenstate that has max mean PR
pr_max_idx = np.argmax(pr_mean, axis=1)
pr_max_sem = pr_sem[np.arange(len(sigmas)), pr_max_idx]

# ── Figure ──
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12, 5))

# ── Panel (a): Clean eigenstate PRs ──
# Color by character: pairs get distinct colors
pair_colors = {
    0: '#1b9e77',  # Trp6-Trp8 pair
    1: '#666666',  # Trp1 isolated
    2: '#666666',  # Trp5 isolated
    3: '#1b9e77',  # Trp6-Trp8 pair
    4: '#d95f02',  # Trp4-Trp7 pair (strongest)
    5: '#d95f02',  # Trp4-Trp7 pair
    6: '#7570b3',  # Trp2-Trp3 pair
    7: '#7570b3',  # Trp2-Trp3 pair
}
colors = [pair_colors[k] for k in range(8)]

bars = ax_a.bar(range(1, 9), pr_clean, color=colors, edgecolor='white', width=0.7)

# Annotate max PR
brightest = np.argmax(pr_clean)
ax_a.text(brightest + 1, pr_clean[brightest] + 0.08,
          f'PR_max = {pr_clean[brightest]:.2f}',
          ha='center', fontsize=12, fontweight='bold', color='#d95f02')

# ax_a.axhline(1.0, color='gray', ls='--', lw=1, alpha=0.5)
# ax_a.axhline(8.0, color='gray', ls=':', lw=1, alpha=0.3)

# ax_a.text(7.5, 1.15, 'PR = 1\n(localised)', fontsize=9, color='gray', ha='right')
# ax_a.text(7.5, 8.15, 'PR = 8\n(fully delocalised)', fontsize=9, color='gray', ha='right')

ax_a.set_xlabel('Eigenstate')
ax_a.set_ylabel(r'Participation ratio PR$_k$')
ax_a.set_title('(a) Clean 8-site Hamiltonian', fontsize=15)
ax_a.set_xticks(range(1, 9))
ax_a.set_xticklabels([r'$\psi_{' + str(i) + r'}$' for i in range(1, 9)])
ax_a.set_ylim(0, 4)

# ── Panel (b): Max PR vs sigma scan ──
ax_b.errorbar(sigmas, pr_max_mean, yerr=pr_max_sem, fmt='o-',
              color='#225ea8', capsize=3, markersize=5, label='Max PR')

# Mark sigma_3
ax_b.axvline(SIGMA_3, color='#d7301f', ls='--', lw=1.5, alpha=0.8)
ax_b.text(SIGMA_3 + 15, 1.8,
          r'$\sigma_3 \approx 200$ cm$^{-1}$' + '\n(MD $T_3$ noise)',
          fontsize=10, color='#d7301f', va='bottom')

# Reference lines
ax_b.axhline(1.0, color='gray', ls='--', lw=1, alpha=0.4)
ax_b.axhline(2.0, color='gray', ls=':', lw=1, alpha=0.3)

ax_b.set_xlabel(r'Added disorder $\sigma$ (cm$^{-1}$)')
ax_b.set_ylabel('Max participation ratio')
ax_b.set_title('(b) Disorder scan (500 realisations)', fontsize=15)
# ax_b.set_ylim(0.8, 2.5)
ax_b.legend(loc='upper right', fontsize=10)

plt.tight_layout()

for ext in ['png', 'pdf']:
    fig.savefig(FIG_DIR / f'fig_participation_ratio.{ext}')
plt.close(fig)
print(f"Saved -> {FIG_DIR}/fig_participation_ratio.{{png,pdf}}")
