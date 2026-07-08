#!/usr/bin/env python
"""Phase 1 / Task 1.2 — full variance decomposition with cross-covariance (RESEARCH_PLAN §2.3).

Linear superposition: ΔE_total = Σ_i ΔE_i  (i ∈ protein, water, nucleotide, ions)
  ⇒ σ_total² = Σ_i σ_i² + 2·Σ_{i<j} Cov(i,j)

Quantifies per-component variance fractions, the protein-water cross-covariance
(focal: dielectric screening), and the aggregate screening ratio
  R_screen = 1 − σ_total² / Σ σ_i²    (R_screen>0 ⇒ destructive interference).

Outputs (results/phase1_basic_stats/):
  variance_decomp.csv     per-site: variance fractions + key cross-cov terms
  screening_ratio.csv     per-site + system-mean screening ratio
  covariance_matrices.npz raw 4×4 covariance + correlation per dataset/site
  variance_decomp.png     stacked σ² + screening gap + 4×4 correlation heatmap
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load, V_TO_CM, TRP_PRETTY, COMPONENTS, phase_dir,
                   write_csv, setup_style, save_fig)


def site_cov(X4):
    """4×4 covariance matrix of the 4 component series for one site.
    X4 : (4, N) array in cm^-1. Returns cov in (cm^-1)²."""
    return np.cov(X4, rowvar=True)   # demeaned by np.cov


def corr_from_cov(C):
    s = np.sqrt(np.diag(C))
    return C / np.outer(s, s)


rows_dec, rows_scr = [], []
cov_store, corr_store = {}, {}
lin_check = []
for tag in ('slow', 'fast'):
    d = load(tag)
    N = d['n_frames']
    # (4, N, 8) cm^-1
    X = np.stack([d[f'delta_s_{c}'] * V_TO_CM for c in COMPONENTS], axis=0)
    sigma_total = d['delta_s_total'].std(axis=0) * V_TO_CM      # (8,)

    cov_mats = np.empty((8, 4, 4))
    corr_mats = np.empty((8, 4, 4))
    for s in range(8):
        C = site_cov(X[:, :, s])
        cov_mats[s] = C
        corr_mats[s] = corr_from_cov(C)

        var_i = np.diag(C)                        # σ_i²  (4,)
        sum_var_i = var_i.sum()
        sum_pair = C[np.triu_indices(4, k=1)].sum() * 2.0   # 2·Σ_{i<j} Cov
        sigma_total_sq = sigma_total[s] ** 2
        sigma_recon_sq = sum_var_i + sum_pair    # = σ_total² by linearity
        lin_check.append((tag, s, sigma_total_sq, sigma_recon_sq))

        screen = 1.0 - sigma_total_sq / sum_var_i if sum_var_i > 0 else float('nan')
        pw_cov = C[0, 1]   # protein-water

        rows_dec.append((
            tag, TRP_PRETTY[s], s,
            f"{var_i[0]:.1f}", f"{var_i[1]:.1f}", f"{var_i[2]:.1f}", f"{var_i[3]:.1f}",
            f"{sum_pair:.1f}",
            *[f"{2*C[i, j] / sigma_total_sq:.4f}" for (i, j) in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]],
            f"{sigma_total_sq:.1f}", f"{sigma_recon_sq:.1f}",
        ))
        rows_scr.append((tag, TRP_PRETTY[s], s,
                         f"{sum_var_i:.1f}", f"{sigma_total_sq:.1f}",
                         f"{screen:.4f}", f"{pw_cov:.1f}",
                         f"{corr_mats[s][0,1]:.4f}"))

    cov_store[tag] = cov_mats
    corr_store[tag] = corr_mats

    # system summary (mean over sites)
    screens = np.array([r[5] for r in rows_scr if r[0] == tag], dtype=float)
    pw_corr = corr_mats[:, 0, 1]
    print(f"[{tag}] N={N}")
    print(f"   σ_total² (cm^-2): mean={sigma_total.mean()**2:.0f}")
    print(f"   screening ratio R_screen: mean={screens.mean():.3f}  "
          f"(min {screens.min():.3f}, max {screens.max():.3f})")
    print(f"   protein-water correlation: mean={pw_corr.mean():.3f}  "
          f"(min {pw_corr.min():.3f}, max {pw_corr.max():.3f})")
    print()

P = phase_dir(1)
# CSV 1: variance decomposition
pair_labels = ['pw', 'pn', 'pi', 'wn', 'wi', 'ni']   # protein-water, protein-nuc, ...
write_csv(P / 'variance_decomp.csv',
          "dataset,site,idx,var_protein,var_water,var_nucleotide,var_ions,"
          "sum_2paircov_cm2,"
          + ",".join(f"frac_2cov_{lab}" for lab in pair_labels) + ","
          "sigma_total_sq_cm2,sigma_recon_sq_cm2",
          rows_dec)
# CSV 2: screening ratio
write_csv(P / 'screening_ratio.csv',
          "dataset,site,idx,sum_var_i_cm2,sigma_total_sq_cm2,screening_ratio,"
          "cov_protein_water_cm2,corr_protein_water",
          rows_scr)
# NPZ: covariance matrices
np.savez_compressed(P / 'covariance_matrices.npz',
                    components=np.array(COMPONENTS),
                    cov_slow=cov_store['slow'], corr_slow=corr_store['slow'],
                    cov_fast=cov_store['fast'], corr_fast=corr_store['fast'])

# linear-superposition check
max_rel_err = max(abs(b - c) / max(b, 1.0) for _, _, b, c in lin_check)
print(f"linearity check max|σ_total² − (Σσ_i²+2ΣCov)|/σ_total² = {max_rel_err:.3e}")
print(f"wrote {P/'variance_decomp.csv'}, screening_ratio.csv, covariance_matrices.npz")

# ------------------------------------------------------------------ figure
setup_style()
import matplotlib.pyplot as plt

# system-averaged correlation matrices (one per dataset)
fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.5),
                        gridspec_kw={'width_ratios': [1, 1, 1.25]})
for ax, tag, title in zip(axs[:2], ('slow', 'fast'), ('slow (50 ns)', 'fast (2 ns)')):
    mean_corr = corr_store[tag].mean(axis=0)
    im = ax.imshow(mean_corr, cmap='RdBu_r', vmin=-1, vmax=1)
    ax.set_xticks(range(4)); ax.set_yticks(range(4))
    short = ['prot', 'water', 'nuc', 'ions']
    ax.set_xticklabels(short, rotation=30, ha='right')
    ax.set_yticklabels(short)
    ax.set_title(f'mean corr  [{title}]\nprot-water r = {mean_corr[0,1]:.2f}')
    for i in range(4):
        for j in range(4):
            ax.text(j, i, f"{mean_corr[i, j]:.2f}", ha='center', va='center',
                    color='white' if abs(mean_corr[i, j]) > 0.55 else 'black',
                    fontsize=8)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

# panel 3: per-site screening gap Σσ_i² vs σ_total²
ax = axs[2]
xs = np.arange(8)
w = 0.38
for off, tag, col in ((-w/2, 'slow', '#225ea8'), (w/2, 'fast', '#41b6c4')):
    sum_var_i = np.array([sum(np.diag(cov_store[tag][s])) for s in range(8)])
    d = load(tag)
    sig_tot_sq = (d['delta_s_total'].std(axis=0) * V_TO_CM) ** 2
    ax.bar(xs + off - w/2, sum_var_i, w*0.9, color=col, alpha=0.45,
           label=f'{tag} · Σσᵢ² (independent sum)')
    ax.bar(xs + off,      sig_tot_sq, w*0.9, color=col,
           label=f'{tag} · σ_total² (actual)')
ax.set_xticks(xs); ax.set_xticklabels(TRP_PRETTY, rotation=20, ha='right')
ax.set_ylabel('variance (cm⁻²)')
ax.set_title('screening gap: Σσᵢ² vs σ_total²\n(bar above α-shaded = cancelled by cross-cov)')
ax.legend(fontsize=8, loc='upper right')

save_fig(fig, P / 'variance_decomp.png')
print(f"wrote {P/'variance_decomp.png'}")
