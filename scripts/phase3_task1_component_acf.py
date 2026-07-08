#!/usr/bin/env python
"""Phase 3 / Task 3.1 — component-wise ACF & timescale fitting (RESEARCH_PLAN §4.2).

Per-component normalized ACF for protein / water / nucleotide / ions, fit to a
single exponential (+ optional plateau) to extract each source's characteristic
correlation time. This is the phase that VALIDATES the Phase 2 τ₁/τ₂/τ₃
assignments: does water own the sub-ps timescales, does protein own the ns one?

Trajectory choice (deviation from plan, documented in PHASE3_REPORT §2):
  * water / nucleotide / ions are fast → fast traj (10 fs) resolves them.
  * protein is slow (ns)               → slow traj (50 ns) resolves it.
  Both trajectories are computed for every component so the choice is visible.

Model:  C(t) = A·exp(−t/τ) + offset      (offset = plateau for non-decayed cases)
τ_int  := A·τ  (integral of the decaying part only; offset is static on this window)

Outputs (results/phase3_source_attribution/):
  comp_acf_fast.npz, comp_acf_slow.npz   raw component ACFs
  comp_tau.csv                           per-component τ (+ A, offset, τ_int)
  comp_acf.png                           4 panels (one per component) fast+slow+fit
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load, V_TO_CM, COMPONENTS, TRP_PRETTY, phase_dir,
                   write_csv, setup_style, save_fig, acf_fft)

T_FAST_MAX_PS = 500.0      # N/4 reliability cutoff (fast)
T_SLOW_MAX_PS = 10_000.0   # N/4 reliability cutoff (slow)
# fit bands: where each traj is trustworthy AND resolves its target dynamics
BAND_FAST = (0.01, 10.0)   # fast-traj band: sub-ps to ~10 ps (above 10 ps slow is better)
BAND_SLOW = (10.0, 8000.0) # slow-traj band: 10 ps to ~8 ns


def exp_plateau(t, A, tau, off):
    """Single exponential decaying from (A+off) to off."""
    return A * np.exp(-t / tau) + off


def biexp(t, A1, A2, t1, t2, off):
    """Bi-exponential + plateau: two decaying components + static offset."""
    return A1 * np.exp(-t / t1) + A2 * np.exp(-t / t2) + off


def fit_one(t, c, band, p0, bounds):
    m = (t >= band[0]) & (t <= band[1])
    if m.sum() < 5:
        return np.full(3, np.nan)
    try:
        p, _ = curve_fit(exp_plateau, t[m], c[m], p0=p0, bounds=bounds,
                         maxfev=40000)
        return p
    except (RuntimeError, ValueError):
        return np.full(3, np.nan)


def main():
    P = phase_dir(3)
    df = load('fast'); ds = load('slow')

    # ---- compute per-component ACFs on both trajectories
    acf_fast = {}   # comp -> (t, acf)
    acf_slow = {}
    for comp in COMPONENTS:
        xf = df[f'delta_s_{comp}'] * V_TO_CM
        xs = ds[f'delta_s_{comp}'] * V_TO_CM
        tf, af = acf_fft(xf, df['dt_ps'])
        tf, af = tf[tf <= T_FAST_MAX_PS], af[tf <= T_FAST_MAX_PS]
        ts, asl = acf_fft(xs, ds['dt_ps'])
        ts, asl = ts[ts <= T_SLOW_MAX_PS], asl[ts <= T_SLOW_MAX_PS]
        acf_fast[comp] = (tf, af)
        acf_slow[comp] = (ts, asl)
    np.savez_compressed(P / 'comp_acf_fast.npz',
                        **{f't_{c}': acf_fast[c][0] for c in COMPONENTS},
                        **{c: acf_fast[c][1] for c in COMPONENTS})
    np.savez_compressed(P / 'comp_acf_slow.npz',
                        **{f't_{c}': acf_slow[c][0] for c in COMPONENTS},
                        **{c: acf_slow[c][1] for c in COMPONENTS})

    # ---- fits: fast band and slow band for each component
    rows = []
    fit_curves = {c: {} for c in COMPONENTS}
    for comp in COMPONENTS:
        tf, af = acf_fast[comp]
        ts, asl = acf_slow[comp]
        # fast-band fit (mean over sites for robustness of the headline τ)
        af_mean = af.mean(axis=1)
        asl_mean = asl.mean(axis=1)
        p_f = fit_one(tf, af_mean, BAND_FAST,
                      p0=[0.5, 0.5, 0.1], bounds=([0, 1e-3, 0], [1.2, 50, 1.0]))
        p_s = fit_one(ts, asl_mean, BAND_SLOW,
                      p0=[0.3, 500.0, 0.05], bounds=([0, 5.0, 0], [1.2, 5e4, 1.0]))
        fit_curves[comp]['fast'] = p_f
        fit_curves[comp]['slow'] = p_s
        A_f, tau_f, off_f = p_f
        A_s, tau_s, off_s = p_s
        rows.append((comp,
                     f"{A_f:.3f}", f"{tau_f:.4f}", f"{off_f:.3f}", f"{A_f*tau_f:.3f}",
                     f"{A_s:.3f}", f"{tau_s:.2f}",  f"{off_s:.3f}", f"{A_s*tau_s:.2f}"))

    write_csv(P / 'comp_tau.csv',
              "component,"
              "A_fast,tau_fast_ps,offset_fast,tau_int_fast_ps,"
              "A_slow,tau_slow_ps,offset_slow,tau_int_slow_ps",
              rows)

    # ---- bonus: bi-exp fast-band fit (validates τ₁ libration vs τ₂ rotation split)
    bi_rows = []
    bi_curves = {}
    for comp in COMPONENTS:
        tf, af = acf_fast[comp]
        m = (tf >= BAND_FAST[0]) & (tf <= BAND_FAST[1])
        y = af.mean(axis=1)[m]
        try:
            p, _ = curve_fit(biexp, tf[m], y,
                             p0=[0.4, 0.3, 0.05, 1.5, 0.2],
                             bounds=([0, 0, 1e-3, 0.1, 0], [1.5, 1.5, 2.0, 50, 1.0]),
                             maxfev=80000)
            A1, A2, t1, t2, off = p
            # sort by τ ascending so τ1 < τ2
            if t1 > t2:
                A1, A2, t1, t2 = A2, A1, t2, t1
        except (RuntimeError, ValueError):
            A1, A2, t1, t2, off = [np.nan] * 5
        bi_curves[comp] = (A1, A2, t1, t2, off)
        bi_rows.append((comp, f"{A1:.3f}", f"{t1:.4f}", f"{A2:.3f}",
                        f"{t2:.3f}", f"{off:.3f}"))
    write_csv(P / 'comp_tau_biexp_fast.csv',
              "component,A1,tau1_ps,A2,tau2_ps,offset", bi_rows)

    print("=== component timescales — single-exp (mean over 8 sites) ===")
    print(f"{'component':<11} {'τ_fast (ps)':>13} {'off_f':>7}   "
          f"{'τ_slow (ps)':>12} {'off_s':>7}")
    for r in rows:
        print(f"{r[0]:<11} {r[2]:>13} {r[3]:>7}   {r[6]:>12} {r[7]:>7}")
    print("\n=== bi-exp fast-band fit (τ₁ libration vs τ₂ rotation split) ===")
    print(f"{'component':<11} {'τ1 (ps)':>9} (A1)   {'τ2 (ps)':>9} (A2)   offset")
    for r in bi_rows:
        print(f"{r[0]:<11} {r[2]:>9} ({r[1]})   {r[4]:>9} ({r[3]})   {r[5]}")
    print("\nPhase 2 reference:  τ₁=0.044 ps (libration)  τ₂=1.70 ps (rotation)  "
          "τ₃=1140–2663 ps (protein)")

    # ---- figure: 4 panels (one per component), fast + slow + fits
    setup_style()
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(10, 6.5), sharex=True)
    colors = {'protein': '#6a51a3', 'water': '#41b6c4',
              'nucleotide': '#fd8d3c', 'ions': '#74c476'}
    for k, comp in enumerate(COMPONENTS):
        ax = axs[k // 2, k % 2]
        tf, af = acf_fast[comp]
        ts, asl = acf_slow[comp]
        ax.semilogx(tf, af.mean(axis=1), color=colors[comp], lw=1.0, alpha=0.55,
                    label=f'fast traj')
        ax.semilogx(ts[ts >= 10], asl[ts >= 10].mean(axis=1),
                    color=colors[comp], lw=1.4, label=f'slow traj')
        # fits
        A_f, tau_f, off_f = fit_curves[comp]['fast']
        A_s, tau_s, off_s = fit_curves[comp]['slow']
        if np.isfinite(tau_f):
            ax.semilogx(tf[(tf >= BAND_FAST[0]) & (tf <= BAND_FAST[1])],
                        exp_plateau(tf[(tf >= BAND_FAST[0]) & (tf <= BAND_FAST[1])],
                                    A_f, tau_f, off_f),
                        'k--', lw=0.8, alpha=0.6)
        if np.isfinite(tau_s):
            ax.semilogx(ts[(ts >= BAND_SLOW[0]) & (ts <= BAND_SLOW[1])],
                        exp_plateau(ts[(ts >= BAND_SLOW[0]) & (ts <= BAND_SLOW[1])],
                                    A_s, tau_s, off_s),
                        'k:', lw=1.0)
        ax.axvline(10.0, color='gray', ls=':', lw=0.6)
        title = f"{comp}"
        if np.isfinite(tau_f) and np.isfinite(tau_s):
            title += f"   τ̂(fast band)={tau_f:.2f} ps   τ̂(slow band)={tau_s:.0f} ps"
        elif np.isfinite(tau_f):
            title += f"   τ̂(fast band)={tau_f:.2f} ps"
        elif np.isfinite(tau_s):
            title += f"   τ̂(slow band)={tau_s:.0f} ps"
        ax.set_title(title, fontsize=8.5)
        ax.set_ylim(-0.1, 1.05)
        if k == 0:
            ax.legend(fontsize=7, loc='upper right')
        if k // 2 == 1:
            ax.set_xlabel('lag t (ps)')
        if k % 2 == 0:
            ax.set_ylabel('C(t)  [self-normalized]')
    fig.suptitle('Phase 3 / Task 3.1 — component-wise ACF (mean over 8 Trp sites)\n'
                 'single-exp fits per band (τ̂ = band-averaged estimate, NOT Phase 2 τ₁/τ₂/τ₃);\n'
                 'dashed = fast-band fit, dotted = slow-band fit; grey : = 10 ps crossover',
                 y=1.0, fontsize=10)
    save_fig(fig, P / 'comp_acf.png')
    print(f"\nwrote comp_acf_fast.npz, comp_acf_slow.npz, comp_tau.csv, comp_acf.png → {P}")


if __name__ == '__main__':
    main()
