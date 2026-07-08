#!/usr/bin/env python
"""Phase 4 / Task 4.1 — 8x8 site cross-correlation matrices (RESEARCH_PLAN §5.2).

Pearson correlation + covariance matrices for total / protein / water /
nucleotide / ions fluctuations, computed on BOTH trajectories.

Trajectory choice rationale:
  * Slow traj (10 ps, 50 ns)  — captures slow conformational fluctuations; this
    is the relevant ensemble for the Phase 5 static-disorder Monte Carlo.
  * Fast traj (10 fs, 2 ns)   — captures fast fluctuations; successive frames
    are autocorrelated so effective-N is smaller, but the instantaneous
    spatial correlation of fast noise is physically informative.

Pair classification (28 unique pairs among 8 sites):
  * intra-A  :  6 pairs among {αW21, αW346, αW388, αW407}  (indices 0-3)
  * intra-B  :  6 pairs among {βW21, βW101, βW344, βW397}  (indices 4-7)
  * inter-AB : 16 pairs crossing the dimer interface

Outputs (results/phase4_spatial_correlation/):
  corr_pearson_slow.npz, corr_pearson_fast.npz    8x8 Pearson r per component
  corr_cov_slow.npz,   corr_cov_fast.npz           8x8 covariance per component
  trp_distances.csv                               pairwise distances (Angstrom)
  corr_summary.csv                                mean |r| by component x class
  corr_heatmap.png                                2x3 heatmaps (slow+fast x total/protein/water)
  corr_vs_distance.png                            |r| vs pair distance
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from itertools import combinations

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load, V_TO_CM, COMPONENTS, TRP_PRETTY, TRP_LABELS, TRP_CHAIN,
                   J_COUPLING, phase_dir, write_csv, setup_style, save_fig)


def pearson_r(X):
    """8x8 Pearson correlation from (N, 8) data matrix. Diagonal=1."""
    Xc = X - X.mean(axis=0)
    cov = (Xc.T @ Xc) / X.shape[0]
    sig = np.sqrt(np.diag(cov))
    r = cov / np.outer(sig, sig)
    np.fill_diagonal(r, 1.0)
    return r, cov, sig


def classify_pair(i, j):
    if TRP_CHAIN[i] == TRP_CHAIN[j]:
        return f'intra-{TRP_CHAIN[i]}'
    return 'inter-AB'


def main():
    P = phase_dir(4)
    df = load('fast'); ds = load('slow')

    # ---- Trp-Trp distances from slow-traj indole positions
    ic = ds['indole_center']
    mean_pos = ic.mean(axis=0)          # (8, 3) Angstrom
    pos_fluc = ic.std(axis=0).mean()    # mean positional RMSF for context
    print(f"mean indole positional RMSF: {pos_fluc:.2f} A")

    pairs_idx = list(combinations(range(8), 2))
    dist = {}
    for i, j in pairs_idx:
        dist[(i, j)] = np.linalg.norm(mean_pos[i] - mean_pos[j])

    write_csv(P / 'trp_distances.csv',
              "site_i,site_j,label_i,label_j,distance_A,classification",
              [(TRP_LABELS[i], TRP_LABELS[j], TRP_PRETTY[i], TRP_PRETTY[j],
                f"{dist[(i,j)]:.2f}", classify_pair(i, j))
               for i, j in pairs_idx])

    nearest_overall  = min(pairs_idx, key=lambda p: dist[p])
    nearest_inter    = min((p for p in pairs_idx if classify_pair(*p) == 'inter-AB'),
                           key=lambda p: dist[p])
    farthest_overall = max(pairs_idx, key=lambda p: dist[p])
    farthest_inter   = max((p for p in pairs_idx if classify_pair(*p) == 'inter-AB'),
                           key=lambda p: dist[p])
    print(f"nearest pair:  {TRP_PRETTY[nearest_overall[0]]}–{TRP_PRETTY[nearest_overall[1]]} "
          f"= {dist[nearest_overall]:.1f} A")
    print(f"nearest inter: {TRP_PRETTY[nearest_inter[0]]}–{TRP_PRETTY[nearest_inter[1]]} "
          f"= {dist[nearest_inter]:.1f} A")
    print(f"farthest:      {TRP_PRETTY[farthest_overall[0]]}–{TRP_PRETTY[farthest_overall[1]]} "
          f"= {dist[farthest_overall]:.1f} A")
    print(f"farthest inter:{TRP_PRETTY[farthest_inter[0]]}–{TRP_PRETTY[farthest_inter[1]]} "
          f"= {dist[farthest_inter]:.1f} A")

    # ---- correlation + covariance matrices on both trajectories
    comps_focus = ['total', 'protein', 'water']
    results = {}
    for tag, d in [('slow', ds), ('fast', df)]:
        for comp in COMPONENTS + ['total']:
            X = d[f'delta_s_{comp}'] * V_TO_CM   # (N, 8) in cm^-1
            r, cov, sig = pearson_r(X)
            results[(tag, comp)] = dict(r=r, cov=cov, sig=sig)

    # save matrices
    for tag in ['slow', 'fast']:
        np.savez_compressed(P / f'corr_pearson_{tag}.npz',
                            **{c: results[(tag, c)]['r'] for c in COMPONENTS + ['total']},
                            labels=np.array(TRP_PRETTY))
        np.savez_compressed(P / f'corr_cov_{tag}.npz',
                            **{c: results[(tag, c)]['cov'] for c in COMPONENTS + ['total']},
                            labels=np.array(TRP_PRETTY),
                            sigma=np.array([results[(tag, c)]['sig'] for c in COMPONENTS + ['total']]))

    # ---- summary: mean |r| by component x classification
    classes = ['intra-A', 'intra-B', 'inter-AB']
    rows = []
    print("\n=== mean |r| by component x class (off-diagonal pairs) ===")
    print(f"{'component':<11} {'traj':>5}  {'intra-A':>8} {'intra-B':>8} {'inter-AB':>9}  {'all':>6}")
    for tag in ['slow', 'fast']:
        for comp in COMPONENTS + ['total']:
            r = results[(tag, comp)]['r']
            by_class = {c: [] for c in classes}
            all_r = []
            for i, j in pairs_idx:
                rij = abs(r[i, j])
                by_class[classify_pair(i, j)].append(rij)
                all_r.append(rij)
            means = {c: np.mean(by_class[c]) for c in classes}
            mean_all = np.mean(all_r)
            rows.append((comp, tag,
                         f"{means['intra-A']:.3f}", f"{means['intra-B']:.3f}",
                         f"{means['inter-AB']:.3f}", f"{mean_all:.3f}"))
            print(f"{comp:<11} {tag:>5}  {means['intra-A']:>8.3f} {means['intra-B']:>8.3f} "
                  f"{means['inter-AB']:>9.3f}  {mean_all:>6.3f}")

    write_csv(P / 'corr_summary.csv',
              "component,trajectory,mean_abs_r_intraA,mean_abs_r_intraB,"
              "mean_abs_r_interAB,mean_abs_r_all", rows)

    # ---- save per-pair table (total/protein/water, slow traj — primary)
    pair_rows = []
    for i, j in pairs_idx:
        d_ij = dist[(i, j)]
        cls = classify_pair(i, j)
        vals = [f"{abs(results[('slow', c)]['r'][i, j]):.3f}" for c in comps_focus]
        pair_rows.append((TRP_PRETTY[i], TRP_PRETTY[j], TRP_LABELS[i], TRP_LABELS[j],
                          f"{d_ij:.1f}", cls, *vals))
    write_csv(P / 'corr_pairs.csv',
              "site_i,site_j,idx_i,idx_j,distance_A,classification,"
              "abs_r_total_slow,abs_r_protein_slow,abs_r_water_slow", pair_rows)

    # ---- condition number of the covariance matrix (for Phase 5 MC)
    print("\n=== covariance matrix properties (Phase 5 MC relevance) ===")
    for tag in ['slow', 'fast']:
        cov = results[(tag, 'total')]['cov']
        eigs = np.linalg.eigvalsh(cov)
        cond = eigs[-1] / eigs[0]
        print(f"total ({tag}): λ_max={eigs[-1]:.0f}  λ_min={eigs[0]:.0f}  "
              f"cond={cond:.1f}  det={np.linalg.det(cov):.2e}")
        print(f"  σ_total (diag sqrt) mean = {np.sqrt(np.diag(cov)).mean():.0f} cm^-1")
    print(f"  reference: J = {J_COUPLING} cm^-1")

    # ---- spatial covariance decomposition: within-component vs cross-component
    # For pair (i,j): Cov(total_i, total_j) = Σ_k Cov(k_i, k_j)
    #                                        + Σ_{k≠l} Cov(k_i, l_j)
    # within-component terms are cooperative (positive); cross-component terms
    # encode dielectric screening (negative) — the spatial analog of Phase 1's
    # screening ratio R_screen.
    print("\n=== spatial covariance decomposition (mean over 28 pairs, slow traj) ===")
    decomp_rows = []
    for tag in ['slow', 'fast']:
        # stack component data: (N, 4, 8)
        Xcomp = np.stack([d[f'delta_s_{c}'] * V_TO_CM
                          for c in COMPONENTS if c != 'total'], axis=-1) \
            if False else None  # unused, see below
        d_ = ds if tag == 'slow' else df
        Xc = {c: d_[f'delta_s_{c}'] * V_TO_CM for c in COMPONENTS}
        within_vals = []; cross_vals = []; total_vals = []
        for i, j in pairs_idx:
            total_ij = 0.0; within_ij = 0.0; cross_ij = 0.0
            for ki in range(len(COMPONENTS)):
                for kj in range(len(COMPONENTS)):
                    ci, cj = COMPONENTS[ki], COMPONENTS[kj]
                    cov_kl = np.mean((Xc[ci][:, i] - Xc[ci][:, i].mean()) *
                                     (Xc[cj][:, j] - Xc[cj][:, j].mean()))
                    total_ij += cov_kl
                    if ki == kj:
                        within_ij += cov_kl
                    else:
                        cross_ij += cov_kl
            within_vals.append(within_ij)
            cross_vals.append(cross_ij)
            total_vals.append(total_ij)
        w, x, t = np.array(within_vals), np.array(cross_vals), np.array(total_vals)
        # report in units of σ_total^2 (mean over sites) for scale
        sig2_mean = results[(tag, 'total')]['sig'].mean() ** 2
        print(f"  {tag}: within = {w.mean():.0f}  cross = {x.mean():.0f}  "
              f"total = {t.mean():.0f}  (cm⁻¹², per pair)")
        print(f"       within/σ² = {w.mean()/sig2_mean:.4f}  "
              f"cross/σ² = {x.mean()/sig2_mean:.4f}  "
              f"total/σ² = {t.mean()/sig2_mean:.4f}")
        print(f"       spatial screening: |cross|/within = {abs(x.mean())/w.mean():.3f}")
        decomp_rows.append((tag, f"{w.mean():.0f}", f"{x.mean():.0f}", f"{t.mean():.0f}",
                            f"{w.mean()/sig2_mean:.4f}", f"{x.mean()/sig2_mean:.4f}",
                            f"{t.mean()/sig2_mean:.4f}",
                            f"{abs(x.mean())/w.mean():.3f}"))
    write_csv(P / 'spatial_cov_decomp.csv',
              "trajectory,within_cm2,cross_cm2,total_cm2,"
              "within_over_sig2,cross_over_sig2,total_over_sig2,"
              "spatial_screening_ratio", decomp_rows)

    # =====================================================================
    # FIGURES
    # =====================================================================
    setup_style()
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    # ---- Figure 0 (headline): single-panel total correlation matrix
    r_total = results[('slow', 'total')]['r']
    fig0, ax0 = plt.subplots(figsize=(6.5, 5.5))
    norm0 = TwoSlopeNorm(vmin=-0.3, vcenter=0, vmax=1.0)
    im0 = ax0.imshow(r_total, cmap='RdBu_r', norm=norm0)
    ax0.set_xticks(range(8)); ax0.set_yticks(range(8))
    ax0.set_xticklabels(TRP_PRETTY, rotation=45, ha='right', fontsize=9)
    ax0.set_yticklabels(TRP_PRETTY, fontsize=9)
    for i in range(8):
        for j in range(8):
            color = 'white' if abs(r_total[i, j]) > 0.4 else 'black'
            ax0.text(j, i, f"{r_total[i, j]:.2f}", ha='center', va='center',
                     fontsize=8, color=color)
    ax0.axhline(3.5, color='k', lw=1.5)
    ax0.axvline(3.5, color='k', lw=1.5)
    fig0.colorbar(im0, ax=ax0, label='Pearson r', shrink=0.8)
    mean_off = np.mean([abs(r_total[i, j]) for i, j in pairs_idx])
    ax0.set_title(f'Total site-energy correlation (slow traj)\n'
                  f'mean off-diagonal |r| = {mean_off:.3f} → negligible',
                  fontsize=10)
    save_fig(fig0, P / 'corr_total.png')

    # ---- Figure 1: 2x3 heatmaps (slow + fast) x (total/protein/water)
    fig, axs = plt.subplots(2, 3, figsize=(13, 9))
    vmin, vmax = -0.5, 1.0
    cmap = 'RdBu_r'
    # diverging norm centred at 0 so sign is visible
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    for row, tag in enumerate(['slow', 'fast']):
        for col, comp in enumerate(comps_focus):
            ax = axs[row, col]
            r = results[(tag, comp)]['r']
            im = ax.imshow(r, cmap=cmap, norm=norm)
            ax.set_xticks(range(8)); ax.set_yticks(range(8))
            ax.set_xticklabels(TRP_PRETTY, rotation=45, ha='right', fontsize=8)
            ax.set_yticklabels(TRP_PRETTY, fontsize=8)
            # annotate cells
            for i in range(8):
                for j in range(8):
                    val = r[i, j]
                    color = 'white' if abs(val) > 0.35 else 'black'
                    ax.text(j, i, f"{val:.2f}", ha='center', va='center',
                            fontsize=6.5, color=color)
            # chain boundary
            ax.axhline(3.5, color='k', lw=1.5, alpha=0.6)
            ax.axvline(3.5, color='k', lw=1.5, alpha=0.6)
            mean_offdiag = np.mean([abs(r[i, j]) for i, j in pairs_idx])
            title = f"{comp} ({tag})\nmean |r| = {mean_offdiag:.3f}"
            ax.set_title(title, fontsize=9)
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.94, 0.15, 0.015, 0.7])
    fig.colorbar(im, cax=cbar_ax, label='Pearson r')
    fig.suptitle('Phase 4 / Task 4.1 — 8×8 site-energy correlation matrices\n'
                 'top: slow traj (50 ns @ 10 ps)   bottom: fast traj (2 ns @ 10 fs)\n'
                 'black lines separate α-chain (0–3) from β-chain (4–7)',
                 y=0.98, fontsize=10)
    save_fig(fig, P / 'corr_heatmap.png')

    # ---- Figure 2: |r| vs distance, slow traj, total/protein/water
    fig2, ax = plt.subplots(figsize=(8, 5.5))
    markers = {'total': 'o', 'protein': 's', 'water': '^'}
    colors_c = {'total': '#333333', 'protein': '#6a51a3', 'water': '#41b6c4'}
    for comp in comps_focus:
        r = results[('slow', comp)]['r']
        xs = [dist[(i, j)] for i, j in pairs_idx]
        ys = [abs(r[i, j]) for i, j in pairs_idx]
        ax.scatter(xs, ys, marker=markers[comp], s=50, alpha=0.7,
                   color=colors_c[comp], label=comp, edgecolors='white', lw=0.5)
    # annotate nearest / farthest
    for p, lbl in [(nearest_overall, 'nearest'), (farthest_overall, 'farthest'),
                   (nearest_inter, 'nearest inter'), (farthest_inter, 'farthest inter')]:
        ax.annotate(f"{lbl}\n{TRP_PRETTY[p[0]]}–{TRP_PRETTY[p[1]]}",
                    xy=(dist[p], abs(results[('slow', 'total')]['r'][p[0], p[1]])),
                    fontsize=6, xytext=(8, 8), textcoords='offset points',
                    color='#555')
    ax.set_xlabel('Trp–Trp distance (Å)')
    ax.set_ylabel('|Pearson r|  (slow traj)')
    ax.set_title('Phase 4 / Task 4.1 — spatial correlation vs distance\n'
                 '(28 unique pairs; slow trajectory)')
    ax.set_ylim(-0.02, 0.6)
    ax.legend(fontsize=9)
    save_fig(fig2, P / 'corr_vs_distance.png')

    # ---- Figure 3: covariance heatmap for total (slow) — what Phase 5 needs
    fig3, axs3 = plt.subplots(1, 2, figsize=(12, 5))
    for ax, comp, title in [(axs3[0], 'total', 'total (slow)'),
                            (axs3[1], 'protein', 'protein (slow)')]:
        cov = results[('slow', comp)]['cov']
        im3 = ax.imshow(cov, cmap='viridis')
        ax.set_xticks(range(8)); ax.set_yticks(range(8))
        ax.set_xticklabels(TRP_PRETTY, rotation=45, ha='right', fontsize=8)
        ax.set_yticklabels(TRP_PRETTY, fontsize=8)
        for i in range(8):
            for j in range(8):
                ax.text(j, i, f"{cov[i,j]:.0f}", ha='center', va='center',
                        fontsize=6.5, color='white' if cov[i, j] > cov.max()*0.5 else 'black')
        ax.axhline(3.5, color='white', lw=1.5, alpha=0.6)
        ax.axvline(3.5, color='white', lw=1.5, alpha=0.6)
        ax.set_title(f"Cov({title})  [cm⁻¹²]\n"
                     f"σ̄ = {np.sqrt(np.diag(cov)).mean():.0f} cm⁻¹", fontsize=9)
        fig3.colorbar(im3, ax=ax, shrink=0.8)
    fig3.suptitle('Phase 4 / Task 4.1 — covariance matrices (Phase 5 MC input)\n'
                  'off-diagonal entries = Cov(δs_i, δs_j) — nonzero ⇒ independent-site '
                  'assumption fails', y=1.0, fontsize=10)
    save_fig(fig3, P / 'cov_heatmap.png')

    print(f"\nwrote corr_heatmap.png, corr_vs_distance.png, cov_heatmap.png, "
          f"corr_summary.csv, corr_pairs.csv, trp_distances.csv → {P}")


if __name__ == '__main__':
    main()
