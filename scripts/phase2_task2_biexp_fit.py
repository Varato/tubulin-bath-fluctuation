#!/usr/bin/env python
"""Phase 2 / Task 2.2 — bi-exponential ACF fitting (RESEARCH_PLAN §3.3).

Model:  C(t) = A_fast·exp(−t/τ_fast) + A_slow·exp(−t/τ_slow),  A_fast+A_slow=1

Two complementary fits per site:
  (1) stitched fit   — bi-exp on the Task 2.1 stitched ACF (plan-compliant headline)
                       + tri-exp comparison for AIC/BIC.
  (2) separate fits  — slow-only single-exp on the slow ACF (t≥10 ps) gives the
                       authoritative τ_slow and A_slow (slow traj captures the
                       true slow-mode variance); fast-only on the fast ACF
                       (t≤10 ps) gives τ_fast.  These interpret the stitched
                       numbers in light of the Task 2.1 sampling-bias finding.

τ_int = A_fast·τ_fast + A_slow·τ_slow  (integrated correlation time).

Outputs (results/phase2_timescale/):
  acf_fit_biexp.csv        per-site stitched-fit params + AIC/BIC (bi vs tri)
  acf_fit_separate.csv     per-site slow-only & fast-only params
  tau_int.csv              per-site τ_int from stitched and separate
  acf_fit.png              stitched ACF + bi-exp fit + components (8 panels)
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load, V_TO_CM, TRP_PRETTY, phase_dir, write_csv,
                   setup_style, save_fig, trapz)


# ---------------- models
def biexp(t, Af, tf, ts):
    return Af * np.exp(-t / tf) + (1 - Af) * np.exp(-t / ts)

def triexp(t, A1, A2, t1, t2, t3):
    return A1 * np.exp(-t / t1) + A2 * np.exp(-t / t2) + (1 - A1 - A2) * np.exp(-t / t3)

def slow_only(t, As, ts):
    return As * np.exp(-t / ts)

def fast_only(t, Af, tf, asymp):
    return Af * np.exp(-t / tf) + asymp


def aic_bic(rss, n, k):
    """AIC and BIC for least-squares fits (Gaussian residuals)."""
    if rss <= 0:
        rss = 1e-30
    aic = n * np.log(rss / n) + 2 * k
    bic = n * np.log(rss / n) + k * np.log(n)
    return aic, bic


def main():
    P = phase_dir(2)
    st = np.load(P / 'acf_stitched.npz')
    t_st, c_st = st['t_lag_ps'], st['acf']          # (N+1, 8), incl t=0
    af = np.load(P / 'acf_fast.npz')
    t_f, c_f = af['t_lag_ps'], af['acf']
    asl = np.load(P / 'acf_slow.npz')
    t_s, c_s = asl['t_lag_ps'], asl['acf']

    # fit masks (exclude t=0 constraint point from stitched fit; incl separately)
    m_st = t_st > 0
    m_f = (t_f > 0) & (t_f <= 10.0)
    m_s = (t_s >= 10.0) & (t_s <= 8000.0)

    bi_rows, sep_rows, ti_rows, tri_rows = [], [], [], []
    fit_curves = {}
    for s in range(8):
        # ---- (1a) bi-exp on stitched
        try:
            p_bi, _ = curve_fit(biexp, t_st[m_st], c_st[m_st, s],
                                p0=[0.7, 0.3, 500.0],
                                bounds=([0.0, 1e-3, 5.0], [1.0, 20.0, 5e4]),
                                maxfev=20000)
            rss_bi = float(((biexp(t_st[m_st], *p_bi) - c_st[m_st, s]) ** 2).sum())
            aic_bi, bic_bi = aic_bic(rss_bi, m_st.sum(), k=3)
        except RuntimeError:
            p_bi = np.array([np.nan] * 3); aic_bi = bic_bi = rss_bi = float('nan')
        Af_bi, tf_bi, ts_bi = p_bi

        # ---- (1b) tri-exp on stitched (AIC/BIC comparison)
        try:
            p_tri, _ = curve_fit(triexp, t_st[m_st], c_st[m_st, s],
                                 p0=[0.4, 0.3, 0.5, 5.0, 500.0],
                                 bounds=([0, 0, 1e-3, 1e-2, 5.0],
                                         [1, 1, 20.0, 1e3, 5e4]),
                                 maxfev=40000)
            rss_tri = float(((triexp(t_st[m_st], *p_tri) - c_st[m_st, s]) ** 2).sum())
            aic_tri, bic_tri = aic_bic(rss_tri, m_st.sum(), k=5)
            # (amp, tau) pairs sorted by tau ascending
            pairs = sorted([(p_tri[0], p_tri[2]),
                            (p_tri[1], p_tri[3]),
                            (1 - p_tri[0] - p_tri[1], p_tri[4])], key=lambda x: x[1])
            amps = np.array([p[0] for p in pairs])
            taus = np.array([p[1] for p in pairs])
        except (RuntimeError, ValueError):
            aic_tri = bic_tri = float('nan')
            amps = np.full(3, np.nan); taus = np.full(3, np.nan)
        dAIC = aic_tri - aic_bi      # >0 favours bi-exp
        dBIC = bic_tri - bic_bi

        # ---- (2a) slow-only fit on slow ACF (authoritative τ_slow, A_slow)
        try:
            p_so, _ = curve_fit(slow_only, t_s[m_s], c_s[m_s, s],
                                p0=[0.3, 1000.0],
                                bounds=([0.0, 50.0], [1.0, 5e4]), maxfev=20000)
        except RuntimeError:
            p_so = np.array([np.nan, np.nan])
        As_slow, ts_slow = p_so

        # ---- (2b) fast-only fit on fast ACF (τ_fast)
        try:
            p_fo, _ = curve_fit(fast_only, t_f[m_f], c_f[m_f, s],
                                p0=[0.6, 0.3, 0.2],
                                bounds=([0.0, 1e-3, 0.0], [1.0, 20.0, 1.0]),
                                maxfev=20000)
        except RuntimeError:
            p_fo = np.array([np.nan] * 3)
        Af_fast, tf_fast, asymp_fast = p_fo

        # τ_int
        ti_bi = Af_bi * tf_bi + (1 - Af_bi) * ts_bi if np.isfinite(Af_bi) else float('nan')
        ti_sep = (1 - As_slow) * tf_fast + As_slow * ts_slow if np.isfinite(As_slow) else float('nan')

        bi_rows.append((TRP_PRETTY[s], s,
                        f"{Af_bi:.4f}", f"{tf_bi:.4f}", f"{ts_bi:.2f}",
                        f"{1 - Af_bi:.4f}", f"{ti_bi:.2f}",
                        f"{aic_bi:.2f}", f"{aic_tri:.2f}", f"{dAIC:+.2f}",
                        f"{bic_bi:.2f}", f"{bic_tri:.2f}", f"{dBIC:+.2f}"))
        sep_rows.append((TRP_PRETTY[s], s,
                         f"{As_slow:.4f}", f"{ts_slow:.2f}",
                         f"{Af_fast:.4f}", f"{tf_fast:.4f}", f"{asymp_fast:.4f}",
                         f"{ti_sep:.2f}"))
        ti_rows.append((TRP_PRETTY[s], s, f"{ti_bi:.2f}", f"{ti_sep:.2f}"))
        tri_rows.append((TRP_PRETTY[s], s,
                         f"{amps[0]:.3f}", f"{taus[0]:.4f}",
                         f"{amps[1]:.3f}", f"{taus[1]:.3f}",
                         f"{amps[2]:.3f}", f"{taus[2]:.1f}"))
        fit_curves[s] = dict(bi=p_bi, slow=p_so, fast=p_fo, tri=p_tri)

    # system means
    def _mean(rows, i):
        v = np.array([float(r[i]) for r in rows if r[i] != 'nan'])
        return v.mean() if len(v) else float('nan')
    print("=== stitched bi-exp fit (system means) ===")
    print(f"   A_fast = {_mean(bi_rows, 2):.3f}   τ_fast = {_mean(bi_rows, 3):.3f} ps")
    print(f"   A_slow = {_mean(bi_rows, 5):.3f}   τ_slow = {_mean(bi_rows, 4):.1f} ps")
    print(f"   τ_int  = {_mean(bi_rows, 6):.2f} ps")
    print(f"   ΔAIC(tri-bi) = {_mean(bi_rows, 9):+.1f}   ΔBIC(tri-bi) = {_mean(bi_rows, 12):+.1f}"
          f"   (positive ⟹ bi-exp preferred)")
    print("\n=== separate fits (system means) ===")
    print(f"   slow-only:  A_slow = {_mean(sep_rows, 2):.3f}   τ_slow = {_mean(sep_rows, 3):.1f} ps")
    print(f"   fast-only:  A_fast = {_mean(sep_rows, 4):.3f}   τ_fast = {_mean(sep_rows, 5):.4f} ps")
    print(f"   τ_int(sep) = {_mean(sep_rows, 7):.2f} ps")

    write_csv(P / 'acf_fit_biexp.csv',
              "site,idx,A_fast,tau_fast_ps,tau_slow_ps,A_slow,tau_int_ps,"
              "AIC_bi,AIC_tri,dAIC_tri_minus_bi,BIC_bi,BIC_tri,dBIC_tri_minus_bi",
              bi_rows)
    write_csv(P / 'acf_fit_separate.csv',
              "site,idx,A_slow_slowfit,tau_slow_slowfit_ps,"
              "A_fast_fastfit,tau_fast_fastfit_ps,asymp_fastfit,tau_int_separate_ps",
              sep_rows)
    write_csv(P / 'tau_int.csv', "site,idx,tau_int_stitched_ps,tau_int_separate_ps", ti_rows)
    write_csv(P / 'acf_fit_triexp.csv',
              "site,idx,A1,tau1_ps,A2,tau2_ps,A3,tau3_ps", tri_rows)
    print("\n=== tri-exp components (sorted by τ, system means) ===")
    print(f"   τ1 = {_mean(tri_rows, 3):.4f} ps  (A1 = {_mean(tri_rows, 2):.3f})   "
          f"[~water libration]")
    print(f"   τ2 = {_mean(tri_rows, 5):.3f} ps  (A2 = {_mean(tri_rows, 4):.3f})   "
          f"[~water rotation]")
    print(f"   τ3 = {_mean(tri_rows, 7):.1f} ps   (A3 = {_mean(tri_rows, 6):.3f})   "
          f"[~protein conformational]")
    print(f"\nwrote acf_fit_biexp.csv, acf_fit_separate.csv, tau_int.csv into {P}")

    # ---- figure
    setup_style()
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 4, figsize=(12, 5.5), sharex=True)
    t_plot = t_st
    for s in range(8):
        ax = axs[s // 4, s % 4]
        ax.semilogx(t_plot[1:], c_st[1:, s], color='k', lw=0.8, alpha=0.5, label='stitched')
        p = fit_curves[s]['bi']
        if np.all(np.isfinite(p)):
            ax.semilogx(t_plot[1:], biexp(t_plot[1:], *p), color='crimson', lw=1.2,
                        label=f'bi-exp  τ_f={p[1]:.2f}ps  τ_s={p[2]:.0f}ps')
            ax.semilogx(t_plot[1:], p[0]*np.exp(-t_plot[1:]/p[1]), color='orange',
                        lw=0.7, ls=':', label=r'$A_f e^{-t/\tau_f}$')
            ax.semilogx(t_plot[1:], (1-p[0])*np.exp(-t_plot[1:]/p[2]), color='purple',
                        lw=0.7, ls=':', label=r'$A_s e^{-t/\tau_s}$')
        ax.axvline(10.0, color='gray', ls=':', lw=0.6)
        ax.set_title(TRP_PRETTY[s], fontsize=9)
        ax.set_ylim(-0.15, 1.05)
        if s == 0:
            ax.legend(fontsize=6.5, loc='upper right')
        if s // 4 == 1:
            ax.set_xlabel('lag t (ps)')
        if s % 4 == 0:
            ax.set_ylabel('C(t)')
    fig.suptitle('Phase 2 / Task 2.2 — bi-exponential fit of stitched ACF', y=1.0)
    save_fig(fig, P / 'acf_fit.png')
    print(f"wrote {P/'acf_fit.png'}")

    # ---- tri-exp figure (three component timescales)
    fig2, axs2 = plt.subplots(2, 4, figsize=(12, 5.5), sharex=True)
    for s in range(8):
        ax = axs2[s // 4, s % 4]
        ax.semilogx(t_plot[1:], c_st[1:, s], color='k', lw=0.8, alpha=0.5,
                    label='stitched ACF')
        p3 = fit_curves[s].get('tri')
        if p3 is not None and np.all(np.isfinite(p3)):
            pairs = sorted([(p3[0], p3[2]), (p3[1], p3[3]),
                            (1 - p3[0] - p3[1], p3[4])], key=lambda x: x[1])
            ax.semilogx(t_plot[1:], triexp(t_plot[1:], *p3), color='crimson',
                        lw=1.2, label='tri-exp total')
            colors = ['#fd8d3c', '#74c476', '#6a51a3']
            labels = [r'$\tau_1$~libration', r'$\tau_2$~rotation', r'$\tau_3$~protein']
            for k, (a, tau) in enumerate(pairs):
                ax.semilogx(t_plot[1:], a * np.exp(-t_plot[1:] / tau),
                            color=colors[k], lw=0.8, ls='--',
                            label=fr'{labels[k]}: $\tau$={tau:.3f} ps, A={a:.2f}')
        ax.axvline(10.0, color='gray', ls=':', lw=0.6)
        ax.set_title(TRP_PRETTY[s], fontsize=9)
        ax.set_ylim(-0.15, 1.05)
        if s == 0:
            ax.legend(fontsize=6, loc='upper right')
        if s // 4 == 1:
            ax.set_xlabel('lag t (ps)')
        if s % 4 == 0:
            ax.set_ylabel('C(t)')
    fig2.suptitle('Phase 2 / Task 2.2 — tri-exponential fit (three physical '
                  'timescales)', y=1.0)
    save_fig(fig2, P / 'acf_fit_triexp.png')
    print(f"wrote {P/'acf_fit_triexp.png'}")


if __name__ == '__main__':
    main()
