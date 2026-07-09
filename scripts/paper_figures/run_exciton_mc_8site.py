#!/usr/bin/env python
"""8-site coloured-noise MC on the full Craddock (2014) Trp network.

Extends the 2-site dimer MC (run_exciton_mc.py) to all eight Trp sites.
Per-site noise is an independent tri-exponential OU sum, justified by the
nearly-diagonal spatial covariance (mean off-diag |r| = 0.019, Appendix A).

  H(t) = H_CRAD + sum_m delta_epsilon_m(t) |m><m|

  delta_epsilon_m(t) = sum_k OU_{m,k}(t),   OU stationary variance = A_k sigma_m^2

Initial excitation on Trp4 (matches lindblad_8site.npz donor=3, the same
donor as the 2-site dimer so the two figures are directly comparable).

Configs: noiseless (H_CRAD only), full (all 3 modes), tau3-only (static
disorder ablation).

Comparison baseline: results/paper_figures/lindblad_8site.npz (Haken-Strobl
pure-dephasing, gamma=50 cm^-1, on the same H_CRAD).

Output: results/paper_figures/exciton_mc_8site.npz
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import CM_TO_RADPS, PAPER_DIR
import run_exciton_mc as M   # reuse H_CRAD, gen_noise, grid, load_site_params

N_SITE = 8
I0 = M.I4   # donor = Trp4 (index 3), matches lindblad_8site.npz


def H_static() -> np.ndarray:
    """Static 8x8 Craddock Hamiltonian in rad/ps (hbar = 1)."""
    return M.H_CRAD * CM_TO_RADPS


def gen_all_noise(sig, tri, rng, keep_modes) -> np.ndarray:
    """Independent tri-exp OU sum per site. Returns (8, N_STEPS) array."""
    out = np.zeros((N_SITE, M.N_STEPS))
    for m in range(N_SITE):
        s = sig[m]
        a, tau = tri[m]
        ka = a.copy()
        for k in range(3):
            if k not in keep_modes:
                ka[k] = 0.0
        out[m] = M.gen_noise(M.N_STEPS, M.DT, s, ka, tau, rng)
    return out


def run_realization(sig, tri, rng, keep_modes) -> np.ndarray:
    """One stochastic Schrödinger realization. Returns P (N_SITE, N_STEPS).

    Exact step propagator (H is frozen within each dt since the noise array is
    sampled at the step points): psi(t+dt) = exp(-i H(t+dt) dt) psi(t).
    """
    noises = gen_all_noise(sig, tri, rng, keep_modes)
    H0 = H_static()
    psi = np.zeros(N_SITE, dtype=complex)
    psi[I0] = 1.0
    P = np.zeros((N_SITE, M.N_STEPS))
    P[:, 0] = np.abs(psi) ** 2
    for i in range(1, M.N_STEPS):
        H = H0.copy()
        H[np.diag_indices(N_SITE)] += noises[:, i]
        psi = expm(-1j * H * M.DT) @ psi
        P[:, i] = np.abs(psi) ** 2
    return P


def run_ensemble(sig, tri, n_real, keep_modes):
    P = np.zeros((n_real, N_SITE, M.N_STEPS))
    t0 = time.time()
    for i in range(n_real):
        rng = np.random.default_rng(42 + i)
        P[i] = run_realization(sig, tri, rng, keep_modes)
        if (i + 1) % 25 == 0 or i == 0:
            rate = (i + 1) / (time.time() - t0)
            print(f"    real {i+1}/{n_real}  ({rate:.2f} Hz, "
                  f"ETA {(n_real - i - 1) / rate / 60:.1f} min)")
    return dict(P=P.mean(0), P_sem=P.std(0) / np.sqrt(n_real))


def run_noiseless() -> np.ndarray:
    """Coherent evolution under H_CRAD only."""
    H0 = H_static()
    psi = np.zeros(N_SITE, dtype=complex)
    psi[I0] = 1.0
    P = np.zeros((N_SITE, M.N_STEPS))
    P[:, 0] = np.abs(psi) ** 2
    U = expm(-1j * H0 * M.DT)
    for i in range(1, M.N_STEPS):
        psi = U @ psi
        P[:, i] = np.abs(psi) ** 2
    return P


def main():
    sig, tri = M.load_site_params()
    print(f"8-site Craddock network. donor = Trp{I0+1} (index {I0}).")
    print(f"Grid: {M.N_STEPS} steps, dt={M.DT*1e3:.0f} fs, T={M.T_MAX} ps.")
    print(f"Per-site sigma (cm^-1): "
          + ", ".join(f"Trp{m+1}={sig[m]:.0f}" for m in range(N_SITE)))
    print()

    out = dict(t=M.TLIST,
               site_names=[f"Trp{m+1}" for m in range(N_SITE)],
               donor=I0, target=M.I7,
               sigma_cm=np.array([sig[m] for m in range(N_SITE)]))

    print("noiseless (H_CRAD only)...")
    P0 = run_noiseless()
    out['P_noiseless'] = P0
    print(f"  P_Trp4(0)={P0[I0,0]:.3f}  P_Trp4(2ps)={P0[I0,-1]:.3f}  "
          f"sum(2ps)={P0[:,-1].sum():.4f}")

    n_real = M.N_REAL   # 500
    for name, keep in [('full', {0, 1, 2}),
                       ('tau3', {2})]:
        print(f"\nMC config '{name}' (modes kept: {sorted(keep)})  "
              f"N_real={n_real}")
        r = run_ensemble(sig, tri, n_real, keep)
        out[f'P_{name}'] = r['P']
        out[f'P_{name}_sem'] = r['P_sem']
        Psum = r['P'].sum(axis=0)
        print(f"  max|sum_m P_m - 1| = {np.max(np.abs(Psum - 1)):.2e}")
        ps = r['P'][:, -1]
        print(f"  P_m(2 ps): "
              + "  ".join(f"Trp{m+1}={ps[m]:.3f}" for m in range(N_SITE)))

    np.savez_compressed(PAPER_DIR / 'exciton_mc_8site.npz', **out)
    print(f"\nsaved -> {PAPER_DIR / 'exciton_mc_8site.npz'}")

    # also compare against existing Lindblad baseline
    lb = np.load(PAPER_DIR / 'lindblad_8site.npz')
    lb_p = lb['P'][:, -1]
    mc_p = out['P_full'][:, -1]
    print("\nP_m(2 ps) comparison (full MC vs Lindblad gamma=50):")
    print(f"  {'site':<6} {'MC':>8} {'Lindblad':>10} {'diff':>8}")
    for m in range(N_SITE):
        print(f"  Trp{m+1:<4} {mc_p[m]:>8.3f} {lb_p[m]:>10.3f} "
              f"{mc_p[m]-lb_p[m]:>+8.3f}")


if __name__ == '__main__':
    main()
