#!/usr/bin/env python
"""Fig 3 (4-panel 2×2): Two-site inhomogeneous dephasing + ENAQT."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import setup, PAPER_DIR, CM_TO_RADPS


def main():
    setup()
    mpl.rcParams['axes.titlesize'] = 17
    mc = np.load(PAPER_DIR / 'exciton_mc.npz')
    t_fs = mc['t'] * 1e3

    fig, axes = plt.subplots(2, 2, figsize=(13, 10), sharey='row')
    plt.subplots_adjust(left=0.08, right=0.96, top=0.92, bottom=0.07,
                        hspace=0.30, wspace=0.08)

    # ── (a) Ablation ──
    ax = axes[0, 0]
    ax.plot(t_fs, mc['P7_noiseless'], ':', color='gray', lw=1.6, alpha=0.7,
            label='noiseless')
    for tag, col, lbl in [
        ('tau3',      '#2CA02C', r'$T_3$ only'),
        ('tau1',      '#d62728', r'$T_1$ only'),
        ('tau1_tau3', '#FF7043', r'$T_1{+}T_3$'),
        ('full',      'black',   'Full')]:
        P7 = mc[f'{tag}_P7']; sem = mc[f'{tag}_P7_sem']
        ax.plot(t_fs, P7, '-', color=col, lw=1.5, label=lbl)
        ax.fill_between(t_fs, P7 - sem, P7 + sem, color=col, alpha=0.12)
    # Lindblad gamma=50 reference
    lb = np.load(PAPER_DIR / 'lindblad_2site.npz')
    ax.plot(lb['t'] * 1e3, lb['P7'], '--', color='#1f77b4', lw=1.7,
            label=r'Lindblad $\gamma{=}50 \text{cm}^{-1}$')
    # White-noise MC verification (single OU: sigma=163, T=5fs -> gamma=50)
    wn = np.load(PAPER_DIR / 'whitenoise_mc.npz')
    ax.plot(wn['t'] * 1e3, wn['P7'], '-', color='#1f77b4', lw=1.5, alpha=0.8,
            label=r'white-noise OU ($\sigma{=}163 \text{cm}^{-1}$, $T{=}5$ fs)')
    ax.set_xlabel('time (fs)'); ax.set_ylabel(r'$P_7(t)$')
    ax.set_title('(a) Noise ablation')
    ax.set_ylim(-0.02, 1.02); ax.set_xlim(0, t_fs[-1])
    ax.legend(loc='upper left', fontsize=9.5, ncol=3)

    # ── (c) T-scan ──
    ax = axes[1, 0]
    enaqt = np.load(PAPER_DIR / 'enaqt_scan.npz')
    tv = enaqt['T_values']; tp = enaqt['T_scan_P7']; ts = enaqt['T_scan_sem']
    p7_t3only = float(mc['tau3_P7'][-1])  # T3-only localization baseline
    ax.errorbar(tv, tp, yerr=ts, fmt='o-', lw=2,
                ms=5, capsize=3)
    ax.axhline(p7_t3only, color='#2CA02C', ls='--', lw=2, alpha=0.6)
    ax.text(tv[-1] * 0.07, p7_t3only + 0.015,
            r'$T_3$ only ($P_7$=' + f'{p7_t3only:.2f})',
            fontsize=9, color='gray')
    # hbar_J = 1.0 / (abs(float(enaqt['J_cm'])) * CM_TO_RADPS)
    # ax.axvline(hbar_J, color='gray', ls='--', lw=2, alpha=0.5)
    # ax.text(hbar_J * 1.3, 0.56, r'$\hbar/J$', fontsize=10, color='gray', fontweight='bold')
    # mark actual MD bath timescales (site-averaged, Table 1)
    md_T = {'$T_1$': (0.046, '#D32F2F'), '$T_2$': (1.77, '#1565C0')}
    for lbl, (Tmd, col) in md_T.items():
        ax.axvline(Tmd, color=col, ls='--', lw=2, alpha=0.8)
        ax.text(Tmd * 0.55, 0.56, lbl, color=col, fontsize=11, fontweight='bold')
    # T3 is off-chart (2663 ps >> 100 ps); annotate at right edge
    # ax.annotate(r'$T_3\to 2663$ ps', xy=(100, 0.56), xytext=(40, 0.56),
    #             fontsize=9, color='#388E3C', fontweight='bold',
    #             arrowprops=dict(arrowstyle='->', color='#388E3C', lw=1.2))
    # ENAQT regime: P7 > 0.3 threshold (well above T3-only ~0.17)
    T_enaqt = 7.2  # ps: where P7 drops to 0.3 on the declining side
    ax.axhline(0.3, color='gray', ls=':', lw=1.2, alpha=0.5)
    ax.axvspan(1e-5, T_enaqt, color='#7B1FA2', alpha=0.06, zorder=0)
    ax.axvline(T_enaqt, color='#7B1FA2', ls='--', lw=1.8, alpha=0.7)
    ax.text(0.0025, 0.31, 'ENAQT', fontsize=10, color='#7B1FA2',
            fontweight='bold')
    ax.text(T_enaqt * 1.2, 0.56, rf'${T_enaqt:.1f}$ ps',
            fontsize=9, color='#7B1FA2', fontweight='bold')
    # Plateau end: T where ENAQT starts declining
    # T_plateau = 0.88  # ps
    # ax.axvline(T_plateau, color='#E65100', ls=':', lw=1.5, alpha=0.7)
    # ax.text(T_plateau * 0.35, 0.40, 'plateau\nend',
    #         fontsize=8, color='#E65100', fontweight='bold')
    ax.set_xscale('log')
    ax.set_xlim(0.5*1e-3, 10000)
    ax.set_xlabel(r'$T_{\mathrm{OU}}$ (ps)')
    ax.set_ylabel(r'$P_7(t = 4\ \mathrm{ps})$')
    ax.set_title(r'(c) Scan correlation time ($\sigma_{\mathrm{OU}}=300$ cm$^{-1}$)')

    # ── (b) Individual traces ──
    ax = axes[0, 1]
    d = np.load(PAPER_DIR / 'mc_realizations.npz')
    P7r = d['P7']
    ax.plot(t_fs, P7r.T, color='gray', alpha=0.3, lw=0.5)
    ax.plot([], [], color='gray', alpha=0.6, lw=1.5, label='individual MC')
    ax.plot(t_fs, mc['full_P7'], color='black', lw=2.2, label='ensemble avg')
    ax.fill_between(t_fs, mc['full_P7'] - mc['full_P7_sem'],
                     mc['full_P7'] + mc['full_P7_sem'], color='black', alpha=0.12)
    ax.set_xlabel('time (fs)')
    # ax.set_ylabel(r'$P_7(t)$')
    ax.set_title('(b) Individual trajectories (Full)')
    # ax.set_ylim(-0.02, 1.0); ax.set_xlim(0, t_fs[-1])
    ax.legend(loc='upper left', fontsize=10, ncol=3)

    # ── (d) sigma-scan ──
    ax = axes[1, 1]
    sv = enaqt['sigma_values']; sp = enaqt['sigma_scan_P7']; ss = enaqt['sigma_scan_sem']
    mask = sv > 0
    ax.errorbar(sv[mask], sp[mask], yerr=ss[mask], fmt='s-',
                lw=2, ms=5, capsize=3)
    ax.axhline(p7_t3only, color='#2CA02C', ls='--', lw=2, alpha=0.6)
    ax.text(sv[mask][-1] * 0.1, p7_t3only + 0.015,
            r'$T_3$ only ($P_7$=' + f'{p7_t3only:.2f})',
            fontsize=9, color='gray')
    # Jcm = float(enaqt['J_cm'])
    # ax.axvline(abs(Jcm), color='gray', ls='--', lw=2, alpha=0.5)
    # ax.text(abs(Jcm) * 1.3, 0.56, r'$J$', fontsize=10, color='gray', fontweight='bold')
    # mark actual MD bath mode amplitudes (site-averaged: sigma * sqrt(f_k))
    md_sigma = {r'$\sigma_1$': (914*np.sqrt(0.51), '#D32F2F', xoff := 1.12),
                r'$\sigma_2$': (914*np.sqrt(0.33), '#1565C0', xoff := 0.66),
                # r'$\sigma_3$': (914*np.sqrt(0.16), '#388E3C')
    }
    for lbl, (smd, col, xoff) in md_sigma.items():
        ax.axvline(smd, color=col, ls='--', lw=2, alpha=0.8)
        ax.text(smd * xoff, 0.56, lbl, color=col, fontsize=11, fontweight='bold')
    # ENAQT window: P7 > 0.3 threshold
    sigma_lo, sigma_hi = 36.0, 9700.0  # cm^-1, where P7 = 0.3
    ax.axhline(0.3, color='gray', ls=':', lw=1.2, alpha=0.5)
    ax.axvspan(sigma_lo, sigma_hi, color='#7B1FA2', alpha=0.06, zorder=0)
    ax.axvline(sigma_lo, color='#7B1FA2', ls='--', lw=1.8, alpha=0.7)
    ax.axvline(sigma_hi, color='#7B1FA2', ls='--', lw=1.8, alpha=0.7)
    ax.text(80, 0.31, 'ENAQT', fontsize=10, color='#7B1FA2',
            fontweight='bold')
    ax.text(sigma_lo * 0.25, 0.56, rf'${sigma_lo:.0f}$ cm$^{{-1}}$',
            fontsize=9, color='#7B1FA2', fontweight='bold')
    ax.text(sigma_hi * 1.2, 0.56, rf'${sigma_hi/1e3:.2f}\times10^3$ cm$^{{-1}}$',
            fontsize=9, color='#7B1FA2', fontweight='bold')
    ax.set_xscale('log')
    ax.set_xlabel(r'$\sigma_{\mathrm{OU}}$ (cm$^{-1}$)')
    # ax.set_ylabel(r'$P_7(4\ \mathrm{ps})$')
    ax.set_title(r'(d) Scan noise amplitude ($T_{\mathrm{OU}}=50$ fs)')
    ax.set_ylim(0, 0.7)

    p_png = PAPER_DIR / 'fig3_exciton.png'
    p_pdf = PAPER_DIR / 'fig3_exciton.pdf'
    fig.savefig(p_png, dpi=300)
    fig.savefig(p_pdf)
    plt.close(fig)
    print(f"saved {p_png}")


if __name__ == '__main__':
    main()
