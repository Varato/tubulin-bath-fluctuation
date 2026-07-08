#!/usr/bin/env python
"""Phase 2 / Task 2.4 — low-sampling resolution bias quantification (RESEARCH_PLAN §3.5).

Downsamples the fast (10 fs) trajectory to 10 ps to mimic a standard MD sampling
rate, then compares τ_int extracted from the downsampled ACF to the true
high-resolution τ_int.  Also fits a bi-exp to the downsampled ACF to show how
τ_fast becomes unresolvable (collapses to ~dt/2).

This is the quantitative justification for the plan's sampling-frequency guideline
(§7.3): sub-ps fluctuations require sub-ps sampling.

Outputs (results/phase2_timescale/):
  sampling_bias.csv        per-site τ_int (true, downsampled), bias factor,
                           τ_fast from downsampled bi-exp fit
  sampling_bias.png        ACF comparison (true vs downsampled) + τ_int bars
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load, V_TO_CM, TRP_PRETTY, phase_dir, write_csv,
                   setup_style, save_fig, acf_fft)

DT_DOWN_PS = 10.0       # target low sampling interval
STRIDE = int(round(DT_DOWN_PS / 0.01))   # 1000


def tau_int_acf(t, c):
    """Integrated correlation time = ∫C(t)dt, truncated at the first zero
    crossing (standard anti-noise truncation; robust to long-lag noise that
    would otherwise contaminate direct integration)."""
    from utils import trapz
    c = np.asarray(c)
    neg = np.where(c <= 0)[0]
    n = int(neg[0]) + 1 if len(neg) else len(c)
    n = max(n, 2)
    return float(trapz(c[:n], t[:n]))


def biexp(t, Af, tf, ts):
    return Af * np.exp(-t / tf) + (1 - Af) * np.exp(-t / ts)


def main():
    P = phase_dir(2)
    df = load('fast')
    x = df['delta_s_total'] * V_TO_CM        # (N, 8) cm⁻¹

    # downsample
    x_ds = x[::STRIDE]                        # (N_ds, 8)
    dt_ds = df['dt_ps'] * STRIDE
    N_ds = x_ds.shape[0]
    print(f"downsampled fast: N={N_ds}  dt={dt_ds:.2f} ps  (span {N_ds*dt_ds/1000:.1f} ns)")

    # ACFs
    t_true, c_true = acf_fft(x, df['dt_ps'])
    t_ds, c_ds = acf_fft(x_ds, dt_ds)
    # truncate at reliability N/4
    t_true = t_true[t_true <= 500.0]; c_true = c_true[:len(t_true)]
    t_ds_max = t_ds[len(t_ds)//4]
    m_ds = t_ds <= t_ds_max
    t_ds, c_ds = t_ds[m_ds], c_ds[m_ds]

    rows = []
    ti_true_arr, ti_ds_arr, tf_ds_arr = np.empty(8), np.empty(8), np.empty(8)
    for s in range(8):
        ti_true = tau_int_acf(t_true, c_true[:, s])
        ti_ds = tau_int_acf(t_ds, c_ds[:, s])
        bias = ti_ds / ti_true if ti_true > 0 else float('nan')
        # bi-exp on downsampled ACF (τ_fast will be unresolvable)
        try:
            p, _ = curve_fit(biexp, t_ds[t_ds > 0], c_ds[t_ds > 0, s],
                             p0=[0.5, DT_DOWN_PS / 2, 500.0],
                             bounds=([0, 1e-3, 5.0], [1.0, DT_DOWN_PS * 2, 5e4]),
                             maxfev=20000)
            tf_ds = float(p[1])
        except RuntimeError:
            tf_ds = float('nan')
        ti_true_arr[s] = ti_true; ti_ds_arr[s] = ti_ds; tf_ds_arr[s] = tf_ds
        rows.append((TRP_PRETTY[s], s, f"{ti_true:.2f}", f"{ti_ds:.2f}",
                     f"{bias:.3f}", f"{tf_ds:.3f}"))

    bias_ratio = ti_ds_arr / ti_true_arr
    print(f"\nτ_int (true, 10 fs):   mean={ti_true_arr.mean():.2f} ps")
    print(f"τ_int (downsamp 10ps): mean={ti_ds_arr.mean():.2f} ps")
    print(f"ratio of means:        {ti_ds_arr.mean()/ti_true_arr.mean():.3f}  "
          f"(cancels fast- and slow-dominated sites)")
    print(f"per-site bias range:   {bias_ratio.min():.2f}–{bias_ratio.max():.2f}  "
          f"(>>1 for fast-dominated sites where the dt/2={DT_DOWN_PS/2:.0f} ps "
          f"\n   trapezoid spike dominates the true sub-ps fast integral)")
    print(f"\nτ_fast from downsampled bi-exp fit: {tf_ds_arr.min():.2f}–"
          f"{tf_ds_arr.max():.2f} ps  (unconstrained — any τ_fast < dt="
          f"{DT_DOWN_PS} ps is fundamentally unresolvable; true τ_fast ≈ "
          f"0.1–0.7 ps from Task 2.2)")
    print(f"\nBottom line: 10 ps sampling destroys all sub-ps timescale info. "
          f"For exciton dephasing (which lives at sub-ps),\n"
          f"this is catastrophic; for τ_int (slow-mode dominated) the mean "
          f"bias is modest but per-site errors are large.")

    write_csv(P / 'sampling_bias.csv',
              "site,idx,tau_int_true_ps,tau_int_downsampled_ps,"
              "bias_factor,tau_fast_downsampled_fit_ps", rows)

    # figure
    setup_style()
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 4, figsize=(12, 5.5))
    for s in range(8):
        ax = axs[s // 4, s % 4]
        ax.semilogx(t_true, c_true[:, s], color='#41b6c4', lw=0.8,
                    label=f'true (10 fs)  τ_int={ti_true_arr[s]:.0f} ps')
        ax.semilogx(t_ds, c_ds[:, s], color='#f03b20', lw=1.0, marker='o',
                    ms=3, label=f'10 ps sampled  τ_int={ti_ds_arr[s]:.0f} ps')
        ax.axhline(0.0, color='gray', ls=':', lw=0.5)
        ax.set_title(TRP_PRETTY[s], fontsize=9)
        ax.set_ylim(-0.2, 1.1)
        ax.set_xlim(5e-3, 2e3)
        if s == 0:
            ax.legend(fontsize=6.5, loc='upper right')
        if s // 4 == 1:
            ax.set_xlabel('lag t (ps)')
        if s % 4 == 0:
            ax.set_ylabel('C(t)')
    fig.suptitle('Phase 2 / Task 2.4 — sampling-rate bias: 10 fs vs 10 ps', y=1.0)
    save_fig(fig, P / 'sampling_bias.png')
    print(f"wrote {P/'sampling_bias.png'}")


if __name__ == '__main__':
    main()
