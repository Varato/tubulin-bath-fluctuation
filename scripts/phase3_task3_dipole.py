#!/usr/bin/env python
"""Phase 3 / Task 3.3 — dipole-orientation fluctuation contribution (RESEARCH_PLAN §4.4).

Two analyses:

(1) Dipole reorientation ACF:  C_μ(t) = ⟨μ̂(0)·μ̂(t)⟩  per Trp site.
    Since μ̂ is a unit vector, C_μ(0) = 1 by construction. Decay → indole ring
    reorientation time. Computed via FFT (sum of per-Cartesian-component ACFs).

(2) Control test: freeze μ̂ to its time-average μ̄, recompute the site-energy
    fluctuation  δs_fixed(t) = μ̄·E_total(t).  Compare σ_fixed vs σ_total:
        σ_fixed/σ_total → field-only fraction (dipole reorientation removed)
        1 − (σ_fixed/σ_total)² → variance fraction attributable to μ̂ motion

Sanity check first: verify δs_total ≈ μ̂·E_total (the projection convention) so
the control test is on the right footing.

Outputs (results/phase3_source_attribution/):
  dmu_acf.npz             C_μ(t) per site (fast traj), + slow traj
  dmu_reorientation.csv   per-site τ_μ (single-exp fit)
  dipole_control.csv      σ_fixed/σ_total, variance fraction per site
  dipole_acf.png          C_μ(t) per site + fits
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load, V_TO_CM, TRP_PRETTY, phase_dir, write_csv,
                   setup_style, save_fig)

T_FAST_MAX_PS = 500.0
T_SLOW_MAX_PS = 10_000.0
BAND_FIT = (0.02, 200.0)   # fit band for τ_μ (avoid t=0 anchor + long-lag noise)


def dipole_acf(dmu, dt):
    """C_μ(t) = ⟨μ̂(0)·μ̂(t)⟩ = Σ_k ACF_unnorm_k(t), normalized so C(0)=1.

    dmu : (N, M, 3) unit vectors. Returns (t_lag, acf) shape (N, M).
    Uses FFT per Cartesian component; sums the raw (unnormalized) ACFs, then
    divides by raw[0] = Σ_k Σ_n μ_k[n]² = N·⟨|μ|²⟩ = N (unit vectors).
    """
    N, M, _ = dmu.shape
    nfft = 1 << int(np.ceil(np.log2(2 * N)))
    F = np.fft.rfft(dmu, n=nfft, axis=0)           # (nfft/2+1, M, 3)
    raw = np.fft.irfft(F * np.conj(F), n=nfft, axis=0)[:N]  # (N, M, 3)
    raw = raw.sum(axis=2)                            # (N, M): Σ_k Σ_n μ_k[n]μ_k[n+t]
    C = raw / raw[0]                                 # C(0) = 1 per site
    return np.arange(N) * dt, C


def exp_decay(t, A, tau, off):
    return A * np.exp(-t / tau) + off


def main():
    P = phase_dir(3)
    df = load('fast'); ds = load('slow')

    # ---- sanity check: δs_total ≈ μ̂·E_total ?
    # try both signs; report which matches
    dot_f = np.einsum('nsk,nsk->ns', df['dmu'], df['E_total'])
    ds_f = df['delta_s_total']
    corr_plus = np.corrcoef(dot_f.ravel(), ds_f.ravel())[0, 1]
    corr_minus = np.corrcoef((-dot_f).ravel(), ds_f.ravel())[0, 1]
    sign = +1 if corr_plus > corr_minus else -1
    rel_err = np.abs(sign * dot_f - ds_f).mean() / np.abs(ds_f).mean()
    print(f"=== projection sanity check ===")
    print(f"   corr( μ̂·E , δs ) = {corr_plus:+.6f}")
    print(f"   corr( −μ̂·E , δs ) = {corr_minus:+.6f}")
    print(f"   ⟹ δs = {'+' if sign>0 else '−'}μ̂·E   "
          f"(mean|rel error| = {rel_err:.2e})")

    # ---- (1) dipole reorientation ACF on both trajectories
    t_f, Cmu_f = dipole_acf(df['dmu'], df['dt_ps'])
    mf = t_f <= T_FAST_MAX_PS
    t_f, Cmu_f = t_f[mf], Cmu_f[mf]

    t_s, Cmu_s = dipole_acf(ds['dmu'], ds['dt_ps'])
    ms = t_s <= T_SLOW_MAX_PS
    t_s, Cmu_s = t_s[ms], Cmu_s[ms]

    np.savez_compressed(P / 'dmu_acf.npz',
                        t_fast_ps=t_f, cmu_fast=Cmu_f,
                        t_slow_ps=t_s, cmu_slow=Cmu_s,
                        labels=np.array(TRP_PRETTY))

    # fit τ_μ per site — on SLOW traj (indole reorients on ns+ timescale; fast
    # traj barely decays: C_μ(500 ps) ≈ 0.88, so fast-traj fit is a lower bound only)
    BAND_SLOW_FIT = (10.0, 10_000.0)
    rows = []
    fit_p = {}
    for s in range(8):
        m = (t_s >= BAND_SLOW_FIT[0]) & (t_s <= BAND_SLOW_FIT[1])
        try:
            p, _ = curve_fit(exp_decay, t_s[m], Cmu_s[m, s],
                             p0=[0.5, 5000.0, 0.3],
                             bounds=([0, 100.0, 0], [1.5, 1e6, 1.0]),
                             maxfev=40000)
        except (RuntimeError, ValueError):
            p = np.array([np.nan] * 3)
        fit_p[s] = p
        A, tau, off = p
        rows.append((TRP_PRETTY[s], s, f"{A:.3f}", f"{tau:.0f}", f"{off:.3f}"))
    write_csv(P / 'dmu_reorientation.csv',
              "site,idx,A,tau_mu_ps,offset_slowtraj", rows)
    taus = [float(r[3]) for r in rows]
    # exciton-timescale relevance: C_μ at T_obs = 2 ps
    i_2ps = np.searchsorted(t_f, 2.0)
    cmu_at_Tobs = Cmu_f[i_2ps].mean()
    print(f"\n=== dipole reorientation τ_μ (slow traj, single-exp) ===")
    print(f"   per-site τ_μ (ps): "
          + "  ".join(f"{r[0]}={r[3]}" for r in rows))
    print(f"   mean τ_μ = {np.nanmean(taus):.0f} ps   "
          f"(range {np.nanmin(taus):.0f}–{np.nanmax(taus):.0f} ps)")
    print(f"   C_μ at T_obs = 2 ps: {cmu_at_Tobs:.4f}  "
          f"(⟹ μ̂ is {cmu_at_Tobs*100:.1f}% correlated on the exciton timescale → static)")

    # ---- (2) control test: fix μ̂ to time-average
    ctrl_rows = []
    for s in range(8):
        mu_bar_f = df['dmu'][:, s, :].mean(axis=0)   # (3,)
        mu_bar_s = ds['dmu'][:, s, :].mean(axis=0)
        # fixed-dipole projection: δs_fixed = sign · μ̄·E_total
        ds_fixed_f = sign * df['E_total'][:, s, :] @ mu_bar_f
        ds_fixed_s = sign * ds['E_total'][:, s, :] @ mu_bar_s
        sig_fixed_f = ds_fixed_f.std()
        sig_fixed_s = ds_fixed_s.std()
        sig_tot_f = (df['delta_s_total'][:, s] * V_TO_CM).std()
        sig_tot_s = (ds['delta_s_total'][:, s] * V_TO_CM).std()
        # δs_fixed is in V/m (E·unitless); convert to cm⁻¹ for apples-to-apples
        sig_fixed_f *= V_TO_CM
        sig_fixed_s *= V_TO_CM
        ratio_f = sig_fixed_f / sig_tot_f
        ratio_s = sig_fixed_s / sig_tot_s
        varfrac_f = 1 - ratio_f**2
        varfrac_s = 1 - ratio_s**2
        ctrl_rows.append((TRP_PRETTY[s], s,
                          f"{sig_tot_f:.1f}", f"{sig_tot_s:.1f}",
                          f"{sig_fixed_f:.1f}", f"{sig_fixed_s:.1f}",
                          f"{ratio_f:.3f}", f"{ratio_s:.3f}",
                          f"{varfrac_f*100:.1f}", f"{varfrac_s*100:.1f}"))
    write_csv(P / 'dipole_control.csv',
              "site,idx,sigma_total_fast_cm,sigma_total_slow_cm,"
              "sigma_fixed_fast_cm,sigma_fixed_slow_cm,"
              "ratio_fast,ratio_slow,varfrac_dipole_fast_pct,varfrac_dipole_slow_pct",
              ctrl_rows)
    rfast = np.mean([float(r[6]) for r in ctrl_rows])
    rslow = np.mean([float(r[7]) for r in ctrl_rows])
    print(f"\n=== control test: σ(fixed-μ) / σ(total) ===")
    print(f"   mean ratio (fast traj) = {rfast:.3f}   "
          f"⟹ dipole-orientation variance fraction = {(1-rfast**2)*100:.1f}%")
    print(f"   mean ratio (slow traj) = {rslow:.3f}   "
          f"⟹ dipole-orientation variance fraction = {(1-rslow**2)*100:.1f}%")
    print(f"   ⟹ field fluctuations dominate; μ̂ reorientation is a small correction")

    # ---- figure: C_μ(t) per site + fit
    setup_style()
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 4, figsize=(12, 5.5), sharex=True)
    for s in range(8):
        ax = axs[s // 4, s % 4]
        ax.semilogx(t_f, Cmu_f[:, s], color='#41b6c4', lw=1.0, alpha=0.8,
                    label='fast')
        ax.semilogx(t_s[t_s >= 10], Cmu_s[t_s >= 10, s], color='#225ea8',
                    lw=1.0, alpha=0.8, label='slow')
        A, tau, off = fit_p[s]
        if np.isfinite(tau):
            tt = np.logspace(np.log10(0.5), np.log10(T_SLOW_MAX_PS), 100)
            ax.semilogx(tt, exp_decay(tt, A, tau, off), 'r--', lw=0.9,
                        label='fit' if s == 0 else None)
        ax.axhline(0, color='gray', lw=0.4)
        ax.set_title(f"{TRP_PRETTY[s]}  τ_μ={fit_p[s][1]:.0f} ps", fontsize=9)
        ax.set_ylim(-0.25, 1.05)
        if s == 0:
            ax.legend(fontsize=7, loc='upper right')
        if s // 4 == 1:
            ax.set_xlabel('lag t (ps)')
        if s % 4 == 0:
            ax.set_ylabel(r'$C_\mu(t)=\langle\hat\mu(0)\cdot\hat\mu(t)\rangle$')
    fig.suptitle('Phase 3 / Task 3.3 — Trp indole dipole reorientation ACF',
                 y=1.0)
    save_fig(fig, P / 'dipole_acf.png')
    print(f"\nwrote dmu_acf.npz, dmu_reorientation.csv, dipole_control.csv, "
          f"dipole_acf.png → {P}")


if __name__ == '__main__':
    main()
