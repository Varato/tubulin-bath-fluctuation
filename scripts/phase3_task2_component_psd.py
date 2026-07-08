#!/usr/bin/env python
"""Phase 3 / Task 3.2 — component-wise PSD spectral decomposition (RESEARCH_PLAN §4.3).

Per-component Welch PSD for protein / water / nucleotide / ions, stitched across
the two trajectories at 1 cm⁻¹ (same recipe as Phase 2 Task 2.3). Then quantify
the band-resolved power-fraction matrix: for each frequency band, what fraction
of the total spectral power comes from each component?

Band definitions (periods shown for intuition; 1 cm⁻¹ ↔ 33.4 ps):
  * slow     : f < 1 cm⁻¹      (periods > 33 ps)   — ns / slow conformational
  * mid      : 1 ≤ f < 50 cm⁻¹  (0.67 – 33 ps)     — ps environmental
  * fast     : 50 ≤ f < 500 cm⁻¹ (0.067 – 0.67 ps) — sub-ps / libration
  * ultrafast: f ≥ 500 cm⁻¹     (< 0.067 ps)        — vibrational tail

NOTE on variance vs PSD: components are cross-correlated (Phase 1: protein-water
r = −0.39), so Σ PSD_component ≠ PSD_total. The fractions below are computed
against Σ_component PSD (the uncorrelated reconstruction), NOT PSD_total. The
difference is the cross-spectral density — reported as a sanity check.

Outputs (results/phase3_source_attribution/):
  comp_psd_fast.npz, comp_psd_slow.npz   per-component Welch PSD
  comp_psd_stitched.npz                  stitched per component
  comp_band_power.csv                     band × component power-fraction matrix
  comp_psd.png                            log-log overlay of 4 components
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load, V_TO_CM, PS_INV_TO_CM, COMPONENTS, TRP_PRETTY,
                   phase_dir, write_csv, setup_style, save_fig, trapz, log_interp)

F_CROSS_CM = 1.0
F_MIN_CM = 0.02
F_MAX_CM = 1668.0
BANDS = [
    ('slow',      0.02,   1.0),
    ('mid',       1.0,   50.0),
    ('fast',     50.0,  500.0),
    ('ultrafast', 500.0, 1668.0),
]


def welch_psd_cm(x, dt_ps, nperseg):
    fs = 1.0 / dt_ps
    f_ps, Pxx = __import__('scipy.signal', fromlist=['welch']).welch(
        x, fs=fs, nperseg=nperseg, detrend='constant', axis=0)
    return f_ps * PS_INV_TO_CM, Pxx / PS_INV_TO_CM


def main():
    P = phase_dir(3)
    df = load('fast'); ds = load('slow')
    nperseg_f = min(df['n_frames'], 16384)
    nperseg_s = min(ds['n_frames'], 2048)

    # ---- per-component PSD on both trajectories
    psd_fast = {}; psd_slow = {}
    for comp in COMPONENTS:
        xf = df[f'delta_s_{comp}'] * V_TO_CM
        xs = ds[f'delta_s_{comp}'] * V_TO_CM
        ff, Sf = welch_psd_cm(xf, df['dt_ps'], nperseg_f)
        fs_, Ss = welch_psd_cm(xs, ds['dt_ps'], nperseg_s)
        psd_fast[comp] = (ff, Sf)
        psd_slow[comp] = (fs_, Ss)
    np.savez_compressed(P / 'comp_psd_fast.npz',
                        f_cm=psd_fast['protein'][0],
                        **{c: psd_fast[c][1] for c in COMPONENTS})
    np.savez_compressed(P / 'comp_psd_slow.npz',
                        f_cm=psd_slow['protein'][0],
                        **{c: psd_slow[c][1] for c in COMPONENTS})

    # ---- stitch per component on common log grid
    f_log = np.logspace(np.log10(F_MIN_CM), np.log10(F_MAX_CM), 300)
    psd_stitch = {}
    for comp in COMPONENTS:
        Sf = log_interp(f_log, *psd_fast[comp])
        Ss = log_interp(f_log, *psd_slow[comp])
        psd_stitch[comp] = np.where((f_log < F_CROSS_CM)[:, None], Ss, Sf)
    np.savez_compressed(P / 'comp_psd_stitched.npz', f_cm=f_log,
                        **{c: psd_stitch[c] for c in COMPONENTS},
                        labels=np.array(TRP_PRETTY))

    # ---- band-resolved power-fraction matrix (mean over sites)
    # per-band, per-component integrated power
    band_comp = np.zeros((len(BANDS), len(COMPONENTS)))
    for bi, (bname, lo, hi) in enumerate(BANDS):
        m = (f_log >= lo) & (f_log <= hi)
        for ci, comp in enumerate(COMPONENTS):
            # mean over 8 sites of the band-integrated power
            band_comp[bi, ci] = trapz(psd_stitch[comp][m].mean(axis=1), f_log[m])
    # fractions: each row (band) sums to 1 over components
    row_sum = band_comp.sum(axis=1, keepdims=True)
    frac = band_comp / row_sum
    # also the absolute power per band (for reporting)
    band_total = band_comp.sum(axis=1)

    rows = []
    for bi, (bname, lo, hi) in enumerate(BANDS):
        rows.append((bname, f"[{lo:.2f},{hi:.0f}] cm⁻¹",
                     f"{band_total[bi]:.2f}",
                     *[f"{frac[bi, ci]*100:.1f}" for ci in range(len(COMPONENTS))]))
    write_csv(P / 'comp_band_power.csv',
              "band,freq_range_cm_inv,total_power_cm2,"
              "protein_pct,water_pct,nucleotide_pct,ions_pct", rows)

    # cross-spectral check: PSD_total vs Σ PSD_component
    xf_tot = df['delta_s_total'] * V_TO_CM
    _, Sf_tot = welch_psd_cm(xf_tot, df['dt_ps'], nperseg_f)
    xs_tot = ds['delta_s_total'] * V_TO_CM
    _, Ss_tot = welch_psd_cm(xs_tot, ds['dt_ps'], nperseg_s)
    S_tot_slow = log_interp(f_log, psd_slow['protein'][0], Ss_tot)
    S_tot_fast = log_interp(f_log, psd_fast['protein'][0], Sf_tot)
    S_tot_stitch = np.where((f_log < F_CROSS_CM)[:, None], S_tot_slow, S_tot_fast)
    S_sum_comp = sum(psd_stitch[c] for c in COMPONENTS)
    # ratio over the full band, mean over sites
    full_mask = np.ones_like(f_log, dtype=bool)
    integ_total = trapz(S_tot_stitch.mean(axis=1), f_log)
    integ_sumcomp = trapz(S_sum_comp.mean(axis=1), f_log)
    print("=== band-resolved component power fractions (mean over 8 sites) ===")
    print(f"{'band':<10} {'range':>18} {'Σpower':>10}  "
          f"protein   water   nucleotide   ions")
    for r in rows:
        print(f"{r[0]:<10} {r[1]:>18} {r[2]:>10}  "
              f"{r[3]:>6}% {r[4]:>6}% {r[5]:>8}% {r[6]:>6}%")
    print(f"\ncross-spectral check (full band):")
    print(f"   ∫PSD_total df      = {integ_total:.1f}  cm⁻¹²")
    print(f"   ∫ΣPSD_comp df      = {integ_sumcomp:.1f}  cm⁻¹²")
    print(f"   ratio Σ_comp/total = {integ_sumcomp/integ_total:.3f}")
    print(f"   (⟨1 if components uncorrelated; the deficit quantifies")
    print(f"    destructive protein-water cross-spectral cancellation⟩)")

    # ---- figure: overlay of 4 component PSDs (mean over sites)
    setup_style()
    import matplotlib.pyplot as plt
    colors = {'protein': '#6a51a3', 'water': '#41b6c4',
              'nucleotide': '#fd8d3c', 'ions': '#74c476'}
    fig, ax = plt.subplots(figsize=(9, 5.5))
    Sw = psd_stitch['water'].mean(axis=1)
    for comp in COMPONENTS:
        ax.loglog(f_log, psd_stitch[comp].mean(axis=1), color=colors[comp],
                  lw=1.3, label=comp)
    ax.loglog(f_log, S_tot_stitch.mean(axis=1), color='k', lw=1.0, alpha=0.5,
              ls='-', label='total (measured)')
    for bname, lo, hi in BANDS:
        ax.axvspan(lo, hi, alpha=0.05, color='gray')
    ax.axvline(F_CROSS_CM, color='gray', ls=':', lw=0.7)

    # annotate decades-below-water at representative frequencies
    for f_rep, label_y_align in [(0.1, 'left'), (10, 'left'), (200, 'left'),
                                  (1000, 'right')]:
        i = np.searchsorted(f_log, f_rep)
        gaps = {}
        for comp in COMPONENTS:
            Sc = psd_stitch[comp].mean(axis=1)
            ratio = Sw[i] / Sc[i] if Sc[i] > 0 else 999
            gaps[comp] = np.log10(ratio) if ratio > 1 else 0
        # build annotation text
        parts = []
        for comp in ['protein', 'ions', 'nucleotide']:
            if gaps[comp] >= 0.1:
                parts.append(f"{comp[:4]}: {gaps[comp]:.1f}d")
        txt = "\n".join(parts)
        y_anchor = Sw[i] * 0.3
        ax.annotate(txt, xy=(f_rep, y_anchor), fontsize=6, color='#555',
                    ha=label_y_align,
                    bbox=dict(boxstyle='round,pad=0.2', fc='white',
                              ec='gray', alpha=0.7))
        # small tick at the reference frequency
        ax.axvline(f_rep, color='gray', ls='-', lw=0.3, alpha=0.3)

    ax.set_xlabel('frequency (cm⁻¹)')
    ax.set_ylabel('PSD (cm⁻¹ / cm⁻¹)  [mean over 8 sites]')
    ax.set_title('Phase 3 / Task 3.2 — component-wise PSD decomposition\n'
                 'annotations = decades below water (d) at representative frequencies',
                 fontsize=9)
    ax.legend(fontsize=8)
    ax.set_ylim(1e-6, None)
    save_fig(fig, P / 'comp_psd.png')

    # ---- per-site breakdown figure (4 bands × 8 sites stacked compactly)
    fig2, axs = plt.subplots(2, 4, figsize=(12, 5.5), sharex=True)
    for s in range(8):
        ax = axs[s // 4, s % 4]
        for comp in COMPONENTS:
            ax.loglog(f_log, psd_stitch[comp][:, s], color=colors[comp],
                      lw=0.8, alpha=0.75, label=comp if s == 0 else None)
        ax.axvline(F_CROSS_CM, color='gray', ls=':', lw=0.6)
        ax.set_title(TRP_PRETTY[s], fontsize=9)
        ax.set_ylim(1e-7, None)
        if s == 0:
            ax.legend(fontsize=6.5, loc='upper right')
        if s // 4 == 1:
            ax.set_xlabel('frequency (cm⁻¹)')
        if s % 4 == 0:
            ax.set_ylabel('PSD (cm⁻¹/cm⁻¹)')
    fig2.suptitle('Phase 3 / Task 3.2 — per-site component PSDs', y=1.0)
    save_fig(fig2, P / 'comp_psd_bysite.png')
    print(f"\nwrote comp_psd.png, comp_psd_bysite.png, comp_band_power.csv → {P}")


if __name__ == '__main__':
    main()
