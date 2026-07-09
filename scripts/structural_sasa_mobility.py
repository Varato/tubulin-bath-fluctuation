#!/usr/bin/env python
"""Merged scatter: τₖ vs whole-residue SASA, colored by dipole mobility θ_loc.

Combines three structural descriptors per Trp site into a single figure:
  - x-axis : whole-residue SASA (solvent exposure)
  - y-axis : τ₁ / τ₂ / τ₃ (relaxation times from Phase 2 tri-exp)
  - color  : θ_loc (local angular deviation from indole_mobility.py)
"""
from __future__ import annotations
import sys, csv, json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import TRP_PRETTY, setup_style, save_fig, RESULTS

OUT = RESULTS / 'structural_correlation'

# ── Whole-residue SASA (nm²) from sasa_analysis/report.md ──
SASA_WHOLE = np.array([0.194, 0.703, 0.046, 0.488, 0.029, 0.066, 0.157, 1.003])

# ── τ from Phase 2 per-site tri-exp ──
def load_tau():
    taus = np.zeros((8, 3))
    with open(OUT.parent / 'phase2_timescale' / 'acf_fit_triexp.csv') as f:
        for row in csv.DictReader(f):
            i = int(row['idx'])
            taus[i] = [float(row['tau1_ps']), float(row['tau2_ps']),
                       float(row['tau3_ps'])]
    return taus

# ── θ_loc from indole_mobility.json ──
def load_mobility():
    with open(OUT / 'indole_mobility.json') as f:
        return np.array(json.load(f)['angular_deviation_local_deg'])


def main():
    setup_style()

    taus = load_tau()
    theta = load_mobility()

    print(f"{'site':<8} {'SASA_whole':>10} {'tau1':>8} {'tau2':>8} {'tau3':>8} "
          f"{'theta_loc':>10}")
    print("-" * 60)
    for i in range(8):
        print(f"{TRP_PRETTY[i]:<8} {SASA_WHOLE[i]:>10.3f} "
              f"{taus[i,0]:>8.4f} {taus[i,1]:>8.3f} {taus[i,2]:>8.1f} "
              f"{theta[i]:>10.1f}")

    tau_labels = [r'$\tau_1$ (libration)',
                  r'$\tau_2$ (rotation)',
                  r'$\tau_3$ (conformational)']

    # ── Figure: 1×3, x=SASA_whole, y=τₖ, color=θ_loc ──
    fig, axs = plt.subplots(1, 3, figsize=(12, 4.5))

    norm = mpl.colors.Normalize(vmin=theta.min(), vmax=theta.max())
    cmap = plt.cm.RdYlBu_r   # red = high mobility, blue = rigid

    for col in range(3):
        ax = axs[col]
        sc = ax.scatter(SASA_WHOLE, taus[:, col], c=theta, cmap=cmap,
                        norm=norm, s=120, zorder=5,
                        edgecolors='black', linewidths=0.6)

        for i in range(8):
            ax.annotate(TRP_PRETTY[i], (SASA_WHOLE[i], taus[i, col]),
                        textcoords='offset points', xytext=(7, 5), fontsize=8)

        # Pearson r
        r = np.corrcoef(SASA_WHOLE, taus[:, col])[0, 1]
        r_log = np.corrcoef(SASA_WHOLE, np.log10(taus[:, col]))[0, 1]
        ax.set_title(f'{tau_labels[col]}\nr = {r:.2f},  r(log) = {r_log:.2f}',
                     fontsize=11)

        ax.set_xlabel('Whole-residue SASA (nm²)', fontsize=11)
        if col == 0:
            ax.set_ylabel('τ (ps)', fontsize=11)
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)

    # Shared colorbar
    cbar = fig.colorbar(sc, ax=axs, orientation='vertical',
                        fraction=0.015, pad=0.02, aspect=40)
    cbar.set_label(r'$\theta_{\mathrm{loc}}$ (deg)', fontsize=11)

    fig.suptitle('Relaxation times vs whole-residue SASA, colored by dipole mobility',
                 fontsize=13, y=1.02)
    save_fig(fig, OUT / 'sasa_tau_mobility_merged.png')
    print(f"\nSaved: {OUT}/sasa_tau_mobility_merged.png")


if __name__ == '__main__':
    main()
