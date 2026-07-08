#!/usr/bin/env python
"""Phase 2 / Task 2.1 — ACF calculation & dual-trajectory stitching (RESEARCH_PLAN §3.2).

Normalized ACF C(t) = <ΔE(0)ΔE(t)> / <ΔE²> via FFT for both trajectories,
stitched into a single full-band ACF.

Stitching rationale (verified empirically):
  * slow traj dt = 10 ps, so it is *unavailable* below 10 ps → fast ACF only
    source for the sub-ps fast-component decay.
  * In the overlap [10 ps, ~500 ps] the two ACFs disagree systematically: the
    2 ns fast trajectory undersamples slow modes, so C_fast(t) decays faster
    than C_slow(t). σ²_slow / σ²_fast ≈ 1.21 (mean). This is the sampling bias
    that Task 2.4 quantifies — slow traj is authoritative for t ≥ 10 ps.
  * Crossover at t = dt_slow = 10 ps. A small discontinuity in the stitched
    C(t) at 10 ps is expected and reported.

Outputs (results/phase2_timescale/):
  acf_fast.npz, acf_slow.npz     raw ACFs (N/2 lags kept)
  acf_stitched.npz               fast(<10 ps) + slow(≥10 ps), log-spaced
  overlap_check.csv              per-site ACF discrepancy in [10,500] ps
  acf_stitched.png               overview
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load, V_TO_CM, TRP_PRETTY, phase_dir, write_csv,
                   setup_style, save_fig, acf_fft, log_interp)

T_FAST_MAX_PS = 500.0     # reliability cutoff N/4 for fast ACF
T_SLOW_MAX_PS = 10_000.0  # reliability cutoff N/4 for slow ACF
T_CROSS_PS = 10.0         # = dt_slow: where slow traj becomes authoritative
OVERLAP = (10.0, 500.0)


def cut_acf(t, acf, t_max):
    m = t <= t_max
    return t[m], acf[m]


def main():
    P = phase_dir(2)
    df = load('fast'); ds = load('slow')
    x_f = df['delta_s_total'] * V_TO_CM
    x_s = ds['delta_s_total'] * V_TO_CM

    t_f, acf_f = acf_fft(x_f, df['dt_ps'])
    t_f, acf_f = cut_acf(t_f, acf_f, T_FAST_MAX_PS)
    np.savez_compressed(P / 'acf_fast.npz', t_lag_ps=t_f, acf=acf_f,
                        dt_ps=df['dt_ps'], labels=np.array(TRP_PRETTY))

    t_s, acf_s = acf_fft(x_s, ds['dt_ps'])
    t_s, acf_s = cut_acf(t_s, acf_s, T_SLOW_MAX_PS)
    np.savez_compressed(P / 'acf_slow.npz', t_lag_ps=t_s, acf=acf_s,
                        dt_ps=ds['dt_ps'], labels=np.array(TRP_PRETTY))

    # variance ratio (the σ² normalization difference between trajectories)
    var_f = (df['delta_s_total'] * V_TO_CM).var(axis=0)
    var_s = (ds['delta_s_total'] * V_TO_CM).var(axis=0)
    var_ratio = var_s / var_f

    # overlap region diagnostic: report mean abs diff + slow/fast ratio of C(t)
    t_ov = np.logspace(np.log10(OVERLAP[0]), np.log10(OVERLAP[1]), 60)
    a_f_ov = log_interp(t_ov, t_f, acf_f)
    a_s_ov = log_interp(t_ov, t_s, acf_s)
    abs_diff = np.abs(a_f_ov - a_s_ov)
    rows = []
    for s in range(8):
        rows.append((TRP_PRETTY[s], s,
                     f"{abs_diff[:, s].mean():.4f}",
                     f"{abs_diff[:, s].max():.4f}",
                     f"{var_ratio[s]:.3f}"))
    write_csv(P / 'overlap_check.csv',
              "site,idx,overlap_mean_abs_diff,overlap_max_abs_diff,var_ratio_slow_over_fast",
              rows)
    print("overlap region [10–500 ps], |C_fast − C_slow|:")
    print(f"   mean over sites: mean={abs_diff.mean():.4f}  max={abs_diff.max():.4f}")
    print(f"   σ²_slow / σ²_fast: mean={var_ratio.mean():.3f}  "
          f"(>1 ⟹ fast traj misses slow-mode variance)")
    print(f"   C_fast/C_slow at 10 ps: "
          f"{(a_f_ov[0].mean()/a_s_ov[0].mean()):.3f}")

    # stitched ACF on log grid
    t_log = np.logspace(np.log10(df['dt_ps']), np.log10(T_SLOW_MAX_PS), 250)
    a_fast_part = log_interp(t_log, t_f, acf_f)
    a_slow_part = log_interp(t_log, t_s, acf_s)
    stitch = np.where((t_log < T_CROSS_PS)[ :, None], a_fast_part, a_slow_part)
    t_st = np.concatenate([[0.0], t_log])
    acf_st = np.vstack([np.ones(8), stitch])
    np.savez_compressed(P / 'acf_stitched.npz', t_lag_ps=t_st, acf=acf_st,
                        t_cross_ps=T_CROSS_PS, labels=np.array(TRP_PRETTY))
    # discontinuity at crossover
    disc = abs(a_fast_part[t_log.searchsorted(T_CROSS_PS)] -
               a_slow_part[t_log.searchsorted(T_CROSS_PS)])
    print(f"stitched ACF discontinuity at {T_CROSS_PS} ps: mean={disc.mean():.4f}")

    # figure
    setup_style()
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 4, figsize=(12, 5.5), sharex=True)
    for s in range(8):
        ax = axs[s // 4, s % 4]
        ax.semilogx(t_f[t_f <= T_FAST_MAX_PS], acf_f[t_f <= T_FAST_MAX_PS, s],
                    color='#41b6c4', lw=0.8, alpha=0.8, label='fast (2 ns)')
        ax.semilogx(t_s[t_s >= OVERLAP[0]], acf_s[t_s >= OVERLAP[0], s],
                    color='#225ea8', lw=0.8, alpha=0.8, label='slow (50 ns)')
        ax.semilogx(t_st[1:], acf_st[1:, s], color='k', lw=1.0, ls='--',
                    label='stitched')
        ax.axvline(T_CROSS_PS, color='gray', ls=':', lw=0.7)
        ax.set_title(TRP_PRETTY[s], fontsize=9)
        ax.set_ylim(-0.15, 1.05)
        if s == 0:
            ax.legend(fontsize=7, loc='upper right')
        if s // 4 == 1:
            ax.set_xlabel('lag t (ps)')
        if s % 4 == 0:
            ax.set_ylabel('C(t)')
    fig.suptitle(f'Phase 2 / Task 2.1 — stitched ACF (crossover {T_CROSS_PS} ps; '
                 f'overlap mean|Δ|={abs_diff.mean():.3f})', y=1.0)
    save_fig(fig, P / 'acf_stitched.png')
    print(f"wrote {P/'acf_stitched.png'}")


if __name__ == '__main__':
    main()
