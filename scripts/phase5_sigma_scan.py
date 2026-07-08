#!/usr/bin/env python
"""Sigma scan: Lindblad vs trajectory MC across weak → strong disorder.

If MC matches Lindblad at small κ and diverges at large κ,
the MC machinery is verified and the divergence is non-Markovian physics.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import qutip as qt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import setup_style, save_fig, phase_dir

CM = 2 * np.pi * 2.998e-2
J_CM = -59.0;  J = J_CM * CM
DT, T_MAX = 0.002, 2.0
TLIST = np.arange(0, T_MAX + DT/2, DT)

k0, k1 = qt.basis(2, 0), qt.basis(2, 1)
N4 = k0*k0.dag(); N7 = k1*k1.dag()
SX = k0*k1.dag() + k1*k0.dag()
H0 = J * SX

A1, TAU1 = 0.50, 0.044
N_REAL = 300

def gen_ou(n, dt, var, tau, rng):
    d = np.exp(-dt / tau)
    ns = np.sqrt(var * (1 - d**2)) * rng.standard_normal(n)
    x = np.empty(n); x[0] = np.sqrt(var) * rng.standard_normal()
    for i in range(1, n): x[i] = x[i-1]*d + ns[i]
    return x

def run_mc(sig_cm, n_real=N_REAL):
    """τ₁-only trajectory MC with equal σ on both sites."""
    sig = sig_cm * CM
    P7 = np.zeros((n_real, len(TLIST)))
    for i in range(n_real):
        rng = np.random.default_rng(i)
        s4 = gen_ou(len(TLIST), DT, A1*sig**2, TAU1, rng)
        s7 = gen_ou(len(TLIST), DT, A1*sig**2, TAU1, rng)
        r = qt.sesolve([H0, [N4, s4], [N7, s7]], k0, TLIST, e_ops=[N7])
        P7[i] = r.expect[0]
    return P7.mean(axis=0)

def run_lindblad(sig_cm):
    """Lindblad mesolve with γ_φ = A₁·2σ²·τ₁ / 2."""
    sig = sig_cm * CM
    g_deph = A1 * 2 * sig**2 * TAU1     # coherence decay rate
    g_L = g_deph / 2                     # collapse op rate
    c_op = np.sqrt(g_L) * qt.sigmaz()
    r = qt.mesolve(H0, k0*k0.dag(), TLIST, c_ops=[c_op], e_ops=[N7])
    return r.expect[0], g_deph


def main():
    P = phase_dir(5); setup_style()

    sigmas = [30, 100, 300, 1000]
    t_fs = TLIST * 1e3

    fig, axs = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    print(f"{'σ(cm⁻¹)':>8} {'κ':>6} {'γ_deph':>8} {'P7_MC':>7} {'P7_LB':>7} {'match?':>8}")
    print("-" * 50)

    for ax, sig_cm in zip(axs.ravel(), sigmas):
        sig = sig_cm * CM
        sig_D = np.sqrt(A1 * 2) * sig           # gap noise std
        kappa = sig_D * TAU1                      # Kubo number

        P7_mc = run_mc(sig_cm)
        P7_lb, g_deph = run_lindblad(sig_cm)

        # noiseless for reference
        r0 = qt.sesolve(H0, k0, TLIST, e_ops=[N7])

        ax.plot(t_fs, r0.expect[0], ':', color='gray', lw=1, alpha=0.3)
        ax.plot(t_fs, P7_lb, '--', color='#377eb8', lw=2, label='Lindblad')
        ax.plot(t_fs, P7_mc, '-', color='#e41a1c', lw=2, label='Trajectory MC')
        ax.fill_between(t_fs,
                        P7_mc - np.zeros_like(P7_mc),  # placeholder
                        P7_mc + np.zeros_like(P7_mc),
                        color='#e41a1c', alpha=0.1)

        match = '✓' if abs(P7_mc[-1] - P7_lb[-1]) < 0.05 else '✗'
        print(f"{sig_cm:8d} {kappa:6.2f} {g_deph:8.1f} "
              f"{P7_mc[-1]:7.3f} {P7_lb[-1]:7.3f} {match:>8}")

        ax.set_title(f'σ = {sig_cm} cm⁻¹,  κ = {kappa:.1f}', fontsize=10)
        ax.set_ylim(-0.02, 1.02)

    for ax in axs[1]: ax.set_xlabel('t (fs)')
    for ax in axs[:, 0]: ax.set_ylabel('P₇(t)')
    axs[0, 0].legend(fontsize=9)
    fig.suptitle('Lindblad vs trajectory MC: agreement at κ≪1, divergence at κ≫1\n'
                 '(gray dotted = noiseless)', y=1.02, fontsize=11)
    save_fig(fig, P / 'sigma_scan_lindblad_vs_mc.png')
    print(f"\nFigure: {P}/sigma_scan_lindblad_vs_mc.png")


if __name__ == '__main__':
    main()
