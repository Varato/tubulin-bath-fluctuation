#!/usr/bin/env python
"""Phase 2 overview figure: stitched ACF (with tri-exp fit) + stitched PSD.

Left:  stitched ACF for all 8 Trps (dots) + tri-exp fit (solid lines)
Right: stitched PSD for all 8 Trps

ACF x-axis: fs.  PSD x-axis: THz (temporal frequency).
"""
from __future__ import annotations
import sys, csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import TRP_PRETTY, setup_style, save_fig, phase_dir

P2 = phase_dir(2)
CM_TO_THZ = 2.998e-2   # 1 cm⁻¹ = 0.02998 THz


def load_triexp():
    """Return (8,3) amps and taus, sorted by τ ascending."""
    amps = np.zeros((8, 3))
    taus = np.zeros((8, 3))
    with open(P2 / 'acf_fit_triexp.csv') as f:
        for row in csv.DictReader(f):
            i = int(row['idx'])
            amps[i] = [float(row['A1']), float(row['A2']), float(row['A3'])]
            taus[i] = [float(row['tau1_ps']), float(row['tau2_ps']),
                       float(row['tau3_ps'])]
    return amps, taus


def main():
    setup_style()

    # ── Load data ──
    acf_data = np.load(P2 / 'acf_stitched.npz')
    t_ps = acf_data['t_lag_ps']
    acf = acf_data['acf']         # (N, 8)

    psd_data = np.load(P2 / 'psd_stitched.npz')
    f_cm = psd_data['f_cm']
    psd = psd_data['psd']         # (N, 8)

    amps, taus = load_triexp()

    # ── Colors ──
    colors = plt.cm.tab10(np.linspace(0, 1, 8))

    t_fs = t_ps * 1e3             # ps → fs
    f_thz = f_cm * CM_TO_THZ      # cm⁻¹ → THz

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # ── Left: ACF + tri-exp fit ──
    for i in range(8):
        # Raw ACF: dots
        ax1.plot(t_fs, acf[:, i], 'o', color=colors[i], ms=2.5,
                 alpha=0.5, markeredgewidth=0)
        # Fit: solid line
        fit = sum(amps[i, k] * np.exp(-t_ps / taus[i, k]) for k in range(3))
        ax1.plot(t_fs, fit, '-', color=colors[i], lw=1.5, alpha=0.9,
                 label=TRP_PRETTY[i])

    ax1.set_xscale('log')
    ax1.set_xlabel('t (fs)')
    ax1.set_ylabel('C(t) / C(0)')
    ax1.set_title('Stitched ACF + tri-exp fit', fontsize=12)
    ax1.set_ylim(-0.05, 1.05)
    ax1.legend(fontsize=8, ncol=2, loc='upper right')

    # ── Right: PSD ──
    for i in range(8):
        ax2.plot(f_thz, psd[:, i], '-', color=colors[i], lw=1.3,
                 alpha=0.85, label=TRP_PRETTY[i])

    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('Frequency (THz)')
    ax2.set_ylabel('PSD')
    ax2.set_title('Stitched PSD', fontsize=12)

    fig.suptitle('Phase 2 — Multi-timescale dynamics overview', fontsize=13,
                 y=1.01)
    plt.tight_layout()
    save_fig(fig, P2 / 'acf_psd_overview.png')
    print(f"Saved: {P2}/acf_psd_overview.png")


if __name__ == '__main__':
    main()
