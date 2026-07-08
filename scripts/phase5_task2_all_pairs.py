#!/usr/bin/env python
"""Phase 5 / Task 5.2 — All-pairs scan with per-site tri-exp parameters.

Uses Craddock 2014 Hamiltonian for J and ΔE, and per-site (A, τ) from
Phase 2 tri-exp fit + per-site σ from Phase 1.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import qutip as qt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import TRP_PRETTY, TRP_LABELS, setup_style, save_fig, phase_dir

CM = 2 * np.pi * 2.998e-2
DT, T_MAX = 0.002, 2.0
TLIST = np.arange(0, T_MAX + DT / 2, DT)
N_STEPS = len(TLIST)
N_REAL = 300
DONOR = 3   # Trp4 = αW407

# ── Craddock 2014 Hamiltonian (cm⁻¹) ──
H_CRAD = np.array([
    [  1,   0, -13,   0,  -2,  -1,   5,  -1],
    [  0, 388, -41,   4,   1,   1,  -4,   1],
    [-13, -41, 342,   2,   0,   1,  -6,   1],
    [  0,   4,   2, 207,  -4,   6, -59,  -1],
    [ -2,   1,   0,  -4,  57,  21,   2,  11],
    [ -1,   1,   1,   6,  21, 102,   5, -51],
    [  5,  -4,  -6, -59,   2,   5, 248,   3],
    [ -1,   1,   1,  -1,  11, -51,   3,   0]
])

# ── Per-site params (Phase 1 σ + Phase 2 per-site tri-exp) ──
#                σ_total  A1     τ1      A2     τ2      A3     τ3
SITE_DATA = [
    ( 843.8, 0.544, 0.0347, 0.376, 0.448,  0.080, 914.3  ),  # 0 αW21
    ( 972.2, 0.453, 0.0257, 0.402, 0.576,  0.145, 459.8  ),  # 1 αW346
    ( 805.8, 0.598, 0.0651, 0.148, 1.165,  0.254, 3882.1 ),  # 2 αW388
    ( 998.5, 0.536, 0.0367, 0.379, 0.698,  0.085, 227.6  ),  # 3 αW407
    ( 720.7, 0.493, 0.0521, 0.123, 0.577,  0.384, 1256.8 ),  # 4 βW21
    ( 846.3, 0.531, 0.0617, 0.317, 8.984,  0.151, 618.1  ),  # 5 βW101
    ( 675.3, 0.508, 0.0326, 0.394, 0.527,  0.098, 308.4  ),  # 6 βW344
    (1453.2, 0.361, 0.0478, 0.524, 0.654,  0.116, 1453.7 ),  # 7 βW397
]
SIGMA   = np.array([s[0] for s in SITE_DATA])
AMPS    = np.array([[s[1], s[3], s[5]] for s in SITE_DATA])   # (8, 3)
TAUS    = np.array([[s[2], s[4], s[6]] for s in SITE_DATA])   # (8, 3) ps

# ablation: which components to keep (1=keep, 0=zero)
CONFIGS = [
    ('full',       [1, 1, 1]),
    ('tau1',       [1, 0, 0]),
    ('tau2',       [0, 1, 0]),
    ('tau3',       [0, 0, 1]),
    ('tau1+tau2',  [1, 1, 0]),
    ('tau1+tau3',  [1, 0, 1]),
]


def gen_ou(n, dt, var, tau, rng):
    d = np.exp(-dt / tau)
    ns = np.sqrt(var * (1 - d**2)) * rng.standard_normal(n)
    x = np.empty(n); x[0] = np.sqrt(var) * rng.standard_normal()
    for i in range(1, n): x[i] = x[i-1]*d + ns[i]
    return x


def gen_site_noise(n, dt, site_idx, keep, rng):
    """Generate noise for a specific site using its per-site (A, τ)."""
    x = np.zeros(n)
    sig = SIGMA[site_idx] * CM
    for k in range(3):
        if keep[k]:
            x += gen_ou(n, dt, AMPS[site_idx, k] * sig**2,
                        TAUS[site_idx, k], rng)
    return x


def run_pair_mc(J_cm, dE_cm, donor, target, keep, n_real):
    """MC for one pair with per-site noise params."""
    J_val = J_cm * CM; dE = dE_cm * CM
    k0, k1 = qt.basis(2, 0), qt.basis(2, 1)
    n0 = k0*k0.dag(); n1 = k1*k1.dag()
    sx = k0*k1.dag() + k1*k0.dag(); sz = n0 - n1
    H0 = J_val * sx + (dE / 2) * sz

    Pa = np.zeros((n_real, N_STEPS))
    for i in range(n_real):
        rng = np.random.default_rng(i)
        s0 = gen_site_noise(N_STEPS, DT, donor, keep, rng)
        s1 = gen_site_noise(N_STEPS, DT, target, keep, rng)
        r = qt.sesolve([H0, [n0, s0], [n1, s1]], k0, TLIST, e_ops=[n1])
        Pa[i] = r.expect[0]
    return Pa.mean(axis=0)


def main():
    P = phase_dir(5); setup_style()
    targets = [i for i in range(8) if i != DONOR]

    print("Per-site parameters:")
    print(f"{'site':<8} {'σ':>6} {'A1':>5} {'τ1':>6} {'A2':>5} {'τ2':>6} "
          f"{'A3':>5} {'τ3':>7}")
    for i in range(8):
        print(f"{TRP_PRETTY[i]:<8} {SIGMA[i]:6.0f} "
              f"{AMPS[i,0]:5.3f} {TAUS[i,0]:6.4f} "
              f"{AMPS[i,1]:5.3f} {TAUS[i,1]:6.3f} "
              f"{AMPS[i,2]:5.3f} {TAUS[i,2]:7.1f}")

    print(f"\nDonor: {TRP_PRETTY[DONOR]} (E={H_CRAD[DONOR,DONOR]})")
    for t in targets:
        J = H_CRAD[DONOR, t]
        dE = H_CRAD[t, t] - H_CRAD[DONOR, DONOR]
        print(f"  → {TRP_PRETTY[t]:<8} J={J:5.0f}  ΔE={dE:5.0f}  "
              f"σ_t={SIGMA[t]:.0f}")

    # ── Run MC (or load cached) ──
    cache = P / 'all_pairs_mc.npz'
    results = {}
    if cache.exists():
        print("Loading cached MC results from all_pairs_mc.npz")
        data = np.load(cache)
        TLIST = data['tlist']
        for cfg_name, _ in CONFIGS:
            for t in targets:
                key = f'{cfg_name}_{TRP_LABELS[t]}'
                if key in data:
                    results[(cfg_name, t)] = data[key]
        print(f"  Loaded {len(results)} curves. Skipping MC.")
    else:
        for cfg_name, keep in CONFIGS:
            print(f"\n  config: {cfg_name}")
            for t in targets:
                J = H_CRAD[DONOR, t]
                dE = H_CRAD[t, t] - H_CRAD[DONOR, DONOR]
                results[(cfg_name, t)] = run_pair_mc(
                    J, dE, DONOR, t, keep, N_REAL)
                print(f"    {TRP_PRETTY[t]:<8} P(2ps)={results[(cfg_name, t)][-1]:.3f}")
        save = dict(tlist=TLIST)
        for (cn, ti), pa in results.items():
            save[f'{cn}_{TRP_LABELS[ti]}'] = pa
        np.savez_compressed(cache, **save)
        print(f"\n  Saved {len(results)} curves → {cache}")

    # ── Save ──
    save = dict(tlist=TLIST)
    for (cn, ti), pa in results.items():
        save[f'{cn}_{TRP_LABELS[ti]}'] = pa
    np.savez_compressed(P / 'all_pairs_mc.npz', **save)

    # ── Plot: single panel, βW344 only (only pair with significant transfer) ──
    t_fs = TLIST * 1e3
    BEST = 6   # βW344
    cfg_colors = {
        'full':       ('black',   '-',  2.0),
        'tau1':       ('#e41a1c', '-',  1.4),
        'tau2':       ('#377eb8', '-',  1.4),
        'tau3':       ('#4daf4a', '-',  1.4),
        'tau1+tau2':  ('#ff7f00', '--', 1.4),
        'tau1+tau3':  ('#984ea3', '--', 1.4),
    }

    fig, ax = plt.subplots(figsize=(8, 5.5))
    # grey line: sum of all other Trps (full config)
    others = sum(results[('full', t)] for t in targets if t != BEST)
    ax.plot(t_fs, others, '-', color='gray', lw=1.5, alpha=0.5,
            label=f'sum of 6 other Trps  (P={others[-1]:.3f})')
    for cfg_name, _ in CONFIGS:
        c, ls, lw = cfg_colors[cfg_name]
        p = results[(cfg_name, BEST)]
        ax.plot(t_fs, p, color=c, ls=ls, lw=lw, alpha=0.85,
                label=f'{cfg_name}  (P={p[-1]:.3f})')

    ax.set_xlabel('t (fs)')
    ax.set_ylabel(r'$P_{\beta\mathrm{W344}}(t)$')
    ax.set_title(
        f'Transfer {TRP_PRETTY[DONOR]} → βW344  '
        f'(J={H_CRAD[DONOR,BEST]:.0f}, ΔE='
        f'{H_CRAD[BEST,BEST]-H_CRAD[DONOR,DONOR]:.0f} cm⁻¹)\n'
        f'All other Trps: P(2ps) < 0.05 (J ≤ 6 cm⁻¹)')
    ax.legend(fontsize=9, loc='lower right')
    ax.set_ylim(-0.02, 0.55)
    save_fig(fig, P / 'all_pairs_population.png')
    print(f"\nFigure: {P}/all_pairs_population.png")


if __name__ == '__main__':
    main()
