#!/usr/bin/env python
"""Compute 8-site MC + Lindblad for all 8 starting Trp sites.

Output: results/paper_figures/fig4_8site_scan.npz
  t:   (N_STEPS,) time in ps
  mc:  (8, 8, N_STEPS) coloured-noise MC populations [start, site, time]
  lb:  (8, 8, N_STEPS) Lindblad populations
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import expm
import qutip as qt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import CM_TO_RADPS, PAPER_DIR
import run_exciton_mc as M

N_SITE = 8
T_MAX = 3.0
DT = M.DT
TLIST = np.arange(0, T_MAX + DT / 2, DT)
N_STEPS = len(TLIST)
N_REAL = 500
GAMMA_CM = 50.0


def run_mc_from_site(sig, tri, start_site, H0, n_real):
    """Colored-noise MC from a given starting site. Returns (8, N_STEPS)."""
    P = np.zeros((n_real, N_SITE, N_STEPS))
    for i in range(n_real):
        rng = np.random.default_rng(42 + i)
        noises = np.zeros((N_SITE, N_STEPS))
        for m in range(N_SITE):
            s = sig[m]
            a, tau = tri[m]
            noises[m] = M.gen_noise(N_STEPS, DT, s, a, tau, rng)
        psi = np.zeros(N_SITE, dtype=complex)
        psi[start_site] = 1.0
        P[i, :, 0] = np.abs(psi) ** 2
        for j in range(1, N_STEPS):
            H = H0.copy()
            H[np.diag_indices(N_SITE)] += noises[:, j] * CM_TO_RADPS
            psi = expm(-1j * H * DT) @ psi
            P[i, :, j] = np.abs(psi) ** 2
    return P.mean(0)


def run_lindblad_from_site(start_site, H_qobj, gamma):
    """Haken–Strobl Lindblad from a given starting site. Returns (8, N_STEPS)."""
    c_ops = [np.sqrt(gamma) * qt.basis(N_SITE, m) * qt.basis(N_SITE, m).dag()
             for m in range(N_SITE)]
    e_ops = [qt.basis(N_SITE, n) * qt.basis(N_SITE, n).dag()
             for n in range(N_SITE)]
    psi0 = qt.basis(N_SITE, start_site)
    r = qt.mesolve(H_qobj, psi0, TLIST, c_ops=c_ops, e_ops=e_ops)
    return np.array(r.expect)


def main():
    sig, tri = M.load_site_params()
    H0 = M.H_CRAD * CM_TO_RADPS
    H_qobj = qt.Qobj(M.H_CRAD * CM_TO_RADPS)
    gamma = GAMMA_CM * CM_TO_RADPS

    print(f"8-site scan: {N_SITE} starting sites, {N_REAL} MC realisations, "
          f"T={T_MAX} ps\n")

    mc_data = {}
    lb_data = {}
    t0 = time.time()
    for s in range(N_SITE):
        print(f"Site Trp{s+1} (index {s})...", flush=True)
        mc_data[s] = run_mc_from_site(sig, tri, s, H0, N_REAL)
        lb_data[s] = run_lindblad_from_site(s, H_qobj, gamma)
        p = mc_data[s]
        print(f"  MC:   P_start={p[s,-1]:.3f}  P_leak={1-p[s,-1]:.3f}")
        p_lb = lb_data[s]
        print(f"  Lind: P_start={p_lb[s,-1]:.3f}  P_leak={1-p_lb[s,-1]:.3f}  "
              f"({time.time()-t0:.0f}s)")

    mc_arr = np.array([mc_data[s] for s in range(N_SITE)])
    lb_arr = np.array([lb_data[s] for s in range(N_SITE)])
    np.savez_compressed(
        PAPER_DIR / 'fig4_8site_scan.npz',
        t=TLIST, mc=mc_arr, lb=lb_arr,
    )
    print(f"\nsaved -> {PAPER_DIR / 'fig4_8site_scan.npz'}")


if __name__ == '__main__':
    main()
