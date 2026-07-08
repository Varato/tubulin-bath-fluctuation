#!/usr/bin/env python
"""Phase 2 / Task 2.3 — full-band PSD calculation & stitching (RESEARCH_PLAN §3.4).

Welch PSD per dataset, converted to cm⁻¹ frequency axis, stitched:
  * f <  1 cm⁻¹  ← slow traj (resolves the slow bath, 0.02–1.67 cm⁻¹)
  * f ≥ 1 cm⁻¹  ← fast traj (resolves the fast bath, up to 1668 cm⁻¹)
Normalization check ∫S df / σ² reported (Phase 0 showed Welch underestimates
by 6–13%; periodogram exact check also reported per site).

Outputs (results/phase2_timescale/):
  psd_fast.npz, psd_slow.npz     Welch PSD (f_cm, psd) per dataset
  psd_stitched.npz               stitched PSD on log-spaced f_cm
  psd_norm_check.csv             per-site ∫Sdf/σ² (stitched Welch + periodogram)
  psd_stitched.png               log-log overview
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.signal import welch, periodogram

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load, V_TO_CM, PS_INV_TO_CM, TRP_PRETTY, phase_dir,
                   write_csv, setup_style, save_fig, trapz, log_interp)

F_CROSS_CM = 1.0      # cm⁻¹, stitch crossover
F_MIN_CM = 0.02       # cm⁻¹, lowest resolvable by slow traj (1/50 ns)
F_MAX_CM = 1668.0     # cm⁻¹, Nyquist of fast traj


def welch_psd_cm(x, dt_ps, nperseg):
    """Welch PSD. Returns (f_cm⁻¹, S_cm) such that ∫S_cm df_cm == σ² up to
    Welch's windowing/demean bias. x: (N,) or (N,M) in cm⁻¹."""
    fs = 1.0 / dt_ps                                   # ps⁻¹
    f_ps, Pxx = welch(x, fs=fs, nperseg=nperseg, detrend='constant', axis=0)
    f_cm = f_ps * PS_INV_TO_CM
    S_cm = Pxx / PS_INV_TO_CM                          # conserve ∫S df
    return f_cm, S_cm


def periodogram_psd_cm(x, dt_ps):
    fs = 1.0 / dt_ps
    f_ps, Pxx = periodogram(x, fs=fs, detrend='constant', window='boxcar', axis=0)
    return f_ps * PS_INV_TO_CM, Pxx / PS_INV_TO_CM


def main():
    P = phase_dir(2)
    df = load('fast'); ds = load('slow')
    x_f = df['delta_s_total'] * V_TO_CM
    x_s = ds['delta_s_total'] * V_TO_CM
    nperseg_f = min(df['n_frames'], 16384)
    nperseg_s = min(ds['n_frames'], 2048)

    f_f, S_f = welch_psd_cm(x_f, df['dt_ps'], nperseg_f)
    f_s, S_s = welch_psd_cm(x_s, ds['dt_ps'], nperseg_s)
    np.savez_compressed(P / 'psd_fast.npz', f_cm=f_f, psd=S_f, labels=np.array(TRP_PRETTY))
    np.savez_compressed(P / 'psd_slow.npz', f_cm=f_s, psd=S_s, labels=np.array(TRP_PRETTY))

    # stitched on log grid
    f_log = np.logspace(np.log10(F_MIN_CM), np.log10(F_MAX_CM), 300)
    Sf_interp = log_interp(f_log, f_f, S_f)
    Ss_interp = log_interp(f_log, f_s, S_s)
    S_stitch = np.where((f_log < F_CROSS_CM)[:, None], Ss_interp, Sf_interp)
    np.savez_compressed(P / 'psd_stitched.npz', f_cm=f_log, psd=S_stitch,
                        f_cross_cm=F_CROSS_CM, labels=np.array(TRP_PRETTY))

    # normalization: stitched Welch + periodogram (exact) per site
    var_f = x_f.var(axis=0); var_s = x_s.var(axis=0)
    rows = []
    for s in range(8):
        # stitched Welch integral (split at crossover for correct f-grid)
        i_slow = f_log < F_CROSS_CM
        i_fast = ~i_slow
        integ_stitch = trapz(S_stitch[i_slow, s], f_log[i_slow]) + \
                       trapz(S_stitch[i_fast, s], f_log[i_fast])
        # periodogram exact on each dataset
        fpo_f, Po_f = periodogram_psd_cm(x_f[:, s], df['dt_ps'])
        fpo_s, Po_s = periodogram_psd_cm(x_s[:, s], ds['dt_ps'])
        peri_fast = trapz(Po_f, fpo_f) / var_f[s]
        peri_slow = trapz(Po_s, fpo_s) / var_s[s]
        rows.append((TRP_PRETTY[s], s,
                     f"{var_f[s] ** 0.5:.2f}", f"{var_s[s] ** 0.5:.2f}",
                     f"{integ_stitch / var_s[s]:.4f}",
                     f"{peri_fast:.5f}", f"{peri_slow:.5f}"))
    print("PSD normalization (∫Sdf/σ²):")
    print(f"   stitched Welch (vs slow σ²): mean="
          f"{np.mean([float(r[4]) for r in rows]):.4f}")
    print(f"   periodogram exact (fast): mean="
          f"{np.mean([float(r[5]) for r in rows]):.5f}")
    print(f"   periodogram exact (slow): mean="
          f"{np.mean([float(r[6]) for r in rows]):.5f}")
    write_csv(P / 'psd_norm_check.csv',
              "site,idx,sigma_fast_cm,sigma_slow_cm,stitched_welch_over_varslow,"
              "periodogram_fast,periodogram_slow", rows)

    # figure
    setup_style()
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 4, figsize=(12, 5.5), sharex=True)
    for s in range(8):
        ax = axs[s // 4, s % 4]
        ax.loglog(f_s, S_s[:, s], color='#225ea8', lw=0.7, alpha=0.7, label='slow')
        ax.loglog(f_f, S_f[:, s], color='#41b6c4', lw=0.7, alpha=0.7, label='fast')
        ax.axvline(F_CROSS_CM, color='gray', ls=':', lw=0.7)
        ax.set_title(TRP_PRETTY[s], fontsize=9)
        ax.set_ylim(1e-8, None)
        if s == 0:
            ax.legend(fontsize=7, loc='upper right')
        if s // 4 == 1:
            ax.set_xlabel('frequency (cm⁻¹)')
        if s % 4 == 0:
            ax.set_ylabel('PSD (cm⁻¹ / cm⁻¹)')
    fig.suptitle(f'Phase 2 / Task 2.3 — stitched PSD (crossover {F_CROSS_CM} cm⁻¹)',
                 y=1.0)
    save_fig(fig, P / 'psd_stitched.png')
    print(f"wrote psd_fast/slow/stitched .npz + psd_norm_check.csv + psd_stitched.png into {P}")


if __name__ == '__main__':
    main()
