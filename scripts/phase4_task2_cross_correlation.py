#!/usr/bin/env python
"""Phase 4 / Task 4.2 — spatial cross-correlation time evolution (§5.3).

Temporal cross-correlation function for representative site pairs, compared
with self-ACF. Answers: how fast do spatial correlations between Trp sites
decay, and does that timescale match the self-ACF?

Cross-correlation definition (normalised):
    C_ij(τ) = ⟨δs_i(t) δs_j(t+τ)⟩ / (σ_i σ_j)

Via FFT (zero-padded to 2N):
    R_ij(τ) = IFFT[ conj(F_i) · F_j ]
which gives Σ_t x_i(t) x_j(t+τ) — correlation of site i "now" with site j "τ later".

Representative pairs (from Task 4.1 distances):
    nearest intra-chain : βW101–βW397  (13.1 Å, intra-B)
    nearest inter-chain : αW407–βW21   (18.0 Å)
    farthest inter-chain: αW346–βW397  (69.4 Å)

Outputs (results/phase4_spatial_correlation/):
    xcorr_total_slow.npz, xcorr_protein_slow.npz   cross-corr for all 28 pairs
    xcorr_water_fast.npz                            water on fast traj (short-lag)
    xcorr_tau.csv                                   τ decay per pair
    xcorr_curves.png                                representative pair curves
    xcorr_fan.png                                   all 28 pairs overlaid
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit
from itertools import combinations

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load, V_TO_CM, COMPONENTS, TRP_PRETTY, TRP_LABELS,
                   phase_dir, write_csv, setup_style, save_fig, log_interp, trapz)

# reliability cutoffs (N/4)
T_FAST_MAX = 500.0    # ps
T_SLOW_MAX = 10_000.0  # ps
T_CROSS_PS = 10.0       # ACF crossover


def xcorr_fft(xi, xj, dt):
    """FFT-based normalised cross-correlation R_ij(τ) = ⟨x_i(t) x_j(t+τ)⟩/(σ_i σ_j).

    Returns (t_lag, xc) where xc has shape (N,) and xc[0] = Pearson r_ij.
    Zero-padded to >= 2N for non-circular correlation.
    """
    xi = np.asarray(xi, float); xj = np.asarray(xj, float)
    xi = xi - xi.mean(); xj = xj - xj.mean()
    N = len(xi)
    nfft = 1 << int(np.ceil(np.log2(2 * N)))
    Fi = np.fft.rfft(xi, n=nfft)
    Fj = np.fft.rfft(xj, n=nfft)
    xc = np.fft.irfft(Fi.conj() * Fj, n=nfft)[:N]
    # normalise by σ_i σ_j N  (so xc[0] = Pearson r)
    si = np.sqrt(np.dot(xi, xi) / N)
    sj = np.sqrt(np.dot(xj, xj) / N)
    xc = xc / (N * si * sj)
    t = np.arange(N) * dt
    return t, xc


def exp_plateau(t, A, tau, off):
    return A * np.exp(-t / tau) + off


def fit_tau(t, c, band, p0, bounds):
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
    P = phase_dir(4)
    df = load('fast'); ds = load('slow')
    pairs = list(combinations(range(8), 2))

    # representative pairs
    rep_pairs = {
        'βW101–βW397\n(13 Å, intra-B)': (5, 7),
        'αW407–βW21\n(18 Å, inter)':    (3, 4),
        'αW346–βW397\n(69 Å, inter)':   (1, 7),
    }

    # ---- compute cross-correlations: total + protein on slow, water on fast
    # Store symmetrized C_sym(τ) = 0.5*[C_ij(τ) + C_ji(τ)] = 0.5*[C_ij(τ) + C_ij(-τ)]
    configs = [
        ('total',   'slow', ds, T_SLOW_MAX),
        ('protein', 'slow', ds, T_SLOW_MAX),
        ('water',   'fast', df, T_FAST_MAX),
    ]
    xcorr_data = {}   # (comp, tag) -> {(i,j): (t, xc_sym)}
    raw_X = {}        # (comp, tag) -> (N, 8) array for on-demand recomputation
    for comp, tag, d, tmax in configs:
        X = d[f'delta_s_{comp}'] * V_TO_CM
        raw_X[(comp, tag)] = (X, d['dt_ps'])
        xc_all = {}
        for i, j in pairs:
            t, xc_ij = xcorr_fft(X[:, i], X[:, j], d['dt_ps'])
            _,  xc_ji = xcorr_fft(X[:, j], X[:, i], d['dt_ps'])
            xc_sym = 0.5 * (xc_ij + xc_ji)
            m = t <= tmax
            xc_all[(i, j)] = (t[m], xc_sym[m])
        xcorr_data[(comp, tag)] = xc_all
        np.savez_compressed(P / f'xcorr_{comp}_{tag}.npz',
                            **{f'{TRP_LABELS[i]}_{TRP_LABELS[j]}_t': xc_all[(i,j)][0]
                               for i, j in pairs},
                            **{f'{TRP_LABELS[i]}_{TRP_LABELS[j]}': xc_all[(i,j)][1]
                               for i, j in pairs})

    # ---- self-ACF for comparison (mean over 8 sites)
    from utils import acf_fft
    self_acf = {}
    for comp, tag, d, tmax in configs:
        X = d[f'delta_s_{comp}'] * V_TO_CM
        dt = d['dt_ps']
        t, ac = acf_fft(X, dt)
        m = t <= tmax
        self_acf[(comp, tag)] = (t[m], ac[m].mean(axis=1))

    # ---- effective sample size (account for autocorrelation)
    # N_eff = N * dt / (2 τ_int)  — for slow traj, N=4001 dt=10ps τ~2500ps → N_eff ~ 8
    # SE(r) ≈ 1/sqrt(N_eff-3).  Pairs below 2*SE are noise.
    print("=== statistical significance of r(0) ===")
    sig_info = {}
    for comp, tag, d, tmax in configs:
        N = d['n_frames']; dt = d['dt_ps']
        # estimate τ_int from self-ACF integral (truncated at first zero)
        t_self, ac_self = self_acf[(comp, tag)]
        # find first zero crossing
        zc = np.where(ac_self < 0)[0]
        tau_int = trapz(ac_self[:zc[0]] if len(zc) > 0 else ac_self, t_self[:zc[0]] if len(zc) > 0 else t_self) \
            if len(zc) > 0 else trapz(ac_self, t_self)
        n_eff = max(N * dt / (2 * tau_int), 3) if tau_int > 0 else N
        se_r = 1.0 / np.sqrt(n_eff - 3)
        sig_info[(comp, tag)] = (n_eff, se_r)
        print(f"  {comp} ({tag}): N={N}, τ_int≈{tau_int:.1f} ps, "
              f"N_eff≈{n_eff:.0f}, SE(r)≈{se_r:.3f}, 2·SE={2*se_r:.3f}")

    # ---- fit decay timescales for all pairs; flag significance
    BAND_SLOW = (10.0, 8000.0)
    BAND_FAST = (0.01, 10.0)
    tau_rows = []
    print("\n=== cross-correlation: r(0) and decay timescale ===")
    print(f"{'pair':<22} {'comp':>9} {'traj':>5} {'r(0)':>8} {'2σ?':>5} {'τ_xc':>10}")
    for comp, tag, d, tmax in configs:
        band = BAND_SLOW if tag == 'slow' else BAND_FAST
        p0 = [0.1, 500.0, 0.0] if tag == 'slow' else [0.05, 0.5, 0.0]
        bounds = ([0, 1.0, -0.5], [1.0, 5e4, 0.5]) if tag == 'slow' \
            else ([0, 1e-3, -0.5], [1.0, 50, 0.5])
        n_eff, se_r = sig_info[(comp, tag)]
        for i, j in pairs:
            t, xc_sym = xcorr_data[(comp, tag)][(i, j)]
            r0 = xc_sym[0]
            significant = abs(r0) > 2 * se_r
            # fit regardless; the fit converges or hits bounds
            p = fit_tau(t, xc_sym, band, p0, bounds)
            A, tau, off = p
            tau_str = f"{tau:.0f}" if np.isfinite(tau) and tau < 5e4 else "—"
            pair_lbl = f"{TRP_PRETTY[i]}–{TRP_PRETTY[j]}"
            tau_rows.append((pair_lbl, TRP_LABELS[i], TRP_LABELS[j], comp, tag,
                             f"{r0:.4f}", f"{n_eff:.0f}", f"{se_r:.3f}",
                             "yes" if significant else "no",
                             tau_str,
                             f"{A:.4f}" if np.isfinite(A) else "nan",
                             f"{off:.4f}" if np.isfinite(off) else "nan"))
            if (i, j) in rep_pairs.values():
                sig_mark = "†" if significant else ""
                print(f"{pair_lbl:<22} {comp:>9} {tag:>5} {r0:>8.4f} {sig_mark:>5} {tau_str:>10}")

    # ---- aggregate significance (mean r over 28 pairs vs SE_mean)
    print("\n=== aggregate significance (mean |r| over 28 pairs) ===")
    for comp, tag, _, _ in configs:
        n_eff, se_r = sig_info[(comp, tag)]
        rs = [abs(xcorr_data[(comp, tag)][(i, j)][1][0]) for i, j in pairs]
        mean_r = np.mean(rs)
        se_mean = se_r / np.sqrt(len(pairs))
        z = mean_r / se_mean
        print(f"  {comp} ({tag}): mean |r| = {mean_r:.4f} ± {se_mean:.4f}  "
              f"→ z = {z:.2f}σ  {'***' if z > 3 else '**' if z > 2 else '*' if z > 1 else 'ns'}")

    # also self-ACF τ for reference
    print("\n=== self-ACF decay (mean over sites, for reference) ===")
    for comp, tag, d, tmax in configs:
        band = BAND_SLOW if tag == 'slow' else BAND_FAST
        p0 = [0.3, 500.0, 0.1] if tag == 'slow' else [0.4, 0.5, 0.2]
        bounds = ([0, 5.0, 0], [1.5, 5e4, 1.0]) if tag == 'slow' \
            else ([0, 1e-3, 0], [1.5, 50, 1.0])
        t, ac = self_acf[(comp, tag)]
        p = fit_tau(t, ac, band, p0, bounds)
        A, tau, off = p
        print(f"  {comp:>9} ({tag:>4}): τ_self = {tau:.1f} ps   A={A:.3f}  off={off:.3f}")

    write_csv(P / 'xcorr_tau.csv',
              "pair,idx_i,idx_j,component,trajectory,r0,N_eff,SE_r,significant,"
              "tau_xc_ps,A,offset", tau_rows)

    # =====================================================================
    # FIGURES
    # =====================================================================
    setup_style()
    import matplotlib.pyplot as plt
    colors_pair = ['#e41a1c', '#377eb8', '#4daf4a']  # red, blue, green
    colors_comp = {'total': '#333333', 'protein': '#6a51a3', 'water': '#41b6c4'}

    # ---- Figure 1: representative pair cross-correlations vs self-ACF
    fig, axs = plt.subplots(1, 3, figsize=(14, 5), sharex=False)
    for col, (comp, tag, _, _) in enumerate(configs):
        ax = axs[col]
        t_self, ac_self = self_acf[(comp, tag)]
        ax.semilogx(t_self, ac_self, color='k', lw=2.0, alpha=0.4,
                    label='self-ACF (mean)')
        for k, (lbl, (i, j)) in enumerate(rep_pairs.items()):
            t, xc_sym = xcorr_data[(comp, tag)][(i, j)]
            pair_short = lbl.split('\n')[0]
            ax.semilogx(t, xc_sym, color=colors_pair[k], lw=1.2,
                        label=pair_short, alpha=0.85)
        if tag == 'slow':
            ax.axvline(T_CROSS_PS, color='gray', ls=':', lw=0.7)
        ax.axhline(0, color='gray', lw=0.5)
        ax.set_ylim(-0.15, 0.5)
        ax.set_title(f"{comp} ({tag} traj)", fontsize=10)
        ax.set_xlabel('lag τ (ps)')
        if col == 0:
            ax.set_ylabel('C_ij(τ)  [symmetrized]')
        ax.legend(fontsize=7, loc='upper right')
    fig.suptitle('Phase 4 / Task 4.2 — spatial cross-correlation vs self-ACF\n'
                 'flat black = self-ACF (site with itself, mean over 8 sites);\n'
                 'coloured = cross-corr for 3 representative pairs',
                 y=1.0, fontsize=9)
    save_fig(fig, P / 'xcorr_curves.png')

    # ---- Figure 2: fan plot — all 28 pairs overlaid (total, slow traj)
    fig2, axs2 = plt.subplots(1, 2, figsize=(12, 5))
    for col, comp in enumerate(['total', 'protein']):
        ax = axs2[col]
        tag = 'slow'
        t_self, ac_self = self_acf[(comp, tag)]
        ax.semilogx(t_self, ac_self, color='k', lw=2.5, alpha=0.3,
                    label='self-ACF')
        for i, j in pairs:
            t, xc_sym = xcorr_data[(comp, tag)][(i, j)]
            ax.semilogx(t, xc_sym, color='steelblue', lw=0.4, alpha=0.4)
        ax.axhline(0, color='gray', lw=0.5)
        ax.axvline(T_CROSS_PS, color='gray', ls=':', lw=0.7)
        ax.set_ylim(-0.15, 0.5)
        ax.set_title(f'{comp} (slow traj)', fontsize=10)
        ax.set_xlabel('lag τ (ps)')
        if col == 0:
            ax.set_ylabel('C_ij(τ)  [symmetrized, 28 pairs]')
        ax.legend(fontsize=8)
    fig2.suptitle('Phase 4 / Task 4.2 — cross-correlation fan plot\n'
                  'all 28 unique pairs overlaid; black = self-ACF',
                 y=1.0, fontsize=9)
    save_fig(fig2, P / 'xcorr_fan.png')

    # ---- Figure 3: r(0) distribution for all pairs (total vs protein, slow traj)
    fig3, ax3 = plt.subplots(figsize=(7, 4.5))
    for comp in ['total', 'protein']:
        r0s = [float(r[5]) for r in tau_rows if r[3] == comp and r[4] == 'slow']
        ax3.hist(r0s, bins=14, alpha=0.5,
                 label=f'{comp} (mean={np.mean(r0s):.3f})',
                 color=colors_comp[comp])
    ax3.axvline(0, color='k', lw=0.8)
    ax3.set_xlabel('r(0) = C_ij(0)  [Pearson correlation]')
    ax3.set_ylabel('count (out of 28 pairs)')
    ax3.set_title('Phase 4 / Task 4.2 — distribution of instantaneous spatial correlation')
    ax3.legend(fontsize=9)
    ax3.set_xlabel('τ_xc (ps)')
    ax3.set_ylabel('count (out of 28 pairs)')
    ax3.set_title('Phase 4 / Task 4.2 — cross-correlation decay timescale distribution')
    ax3.legend(fontsize=9)
    save_fig(fig3, P / 'xcorr_tau_hist.png')

    print(f"\nwrote xcorr_curves.png, xcorr_fan.png, xcorr_tau_hist.png, "
          f"xcorr_tau.csv → {P}")


if __name__ == '__main__':
    main()
