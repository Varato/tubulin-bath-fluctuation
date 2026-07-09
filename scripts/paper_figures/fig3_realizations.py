#!/usr/bin/env python
"""Diagnostic: individual MC realizations of the 2-site dimer under full MD noise.

Runs N_show stochastic Schrödinger realisations of the Trp4-Trp7 dimer driven
by the full tri-exponential OU bath, and plots each P_7(t) trace individually
(thin, translucent) overlaid with the ensemble average (thick black) and the
500-realisation average from exciton_mc.npz (dashed red).

Purpose: inspect whether individual trajectories are coherent (oscillating at
~J) or already smooth, and how the monotonic ensemble average emerges from
dephasing across realisations. This is the diagnostic for whether the transfer
is FRET-like (incoherent, dephasing-washed) at the trajectory level.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import setup, save, CM_TO_RADPS, PAPER_DIR
import run_exciton_mc as M

N_SHOW = 20


def run_one(s4, s7, a4, t4, a7, t7, H0, seed):
    """One stochastic Schrödinger realisation. Returns P7(t), P4(t)."""
    rng = np.random.default_rng(seed)
    n4 = M.gen_noise(M.N_STEPS, M.DT, s4, a4, t4, rng)
    n7 = M.gen_noise(M.N_STEPS, M.DT, s7, a7, t7, rng)
    psi = np.array([1.0, 0.0], dtype=complex)   # start on Trp4
    P4 = np.zeros(M.N_STEPS); P7 = np.zeros(M.N_STEPS)
    P4[0] = abs(psi[0])**2; P7[0] = abs(psi[1])**2
    for i in range(1, M.N_STEPS):
        H = H0 + np.diag([n4[i], n7[i]])
        psi = expm(-1j * H * M.DT) @ psi
        P4[i] = abs(psi[0])**2; P7[i] = abs(psi[1])**2
    return P4, P7


def main():
    setup()
    sig, tri = M.load_site_params()
    s4, s7 = sig[M.I4], sig[M.I7]
    a4, t4 = tri[M.I4]; a7, t7 = tri[M.I7]
    J  = M.J_CM   * CM_TO_RADPS
    dE = M.D_E_CM * CM_TO_RADPS
    H0 = J * np.array([[0, 1], [1, 0]], dtype=complex) \
       + (dE / 2) * np.array([[1, 0], [0, -1]], dtype=complex)

    P7_traces = np.zeros((N_SHOW, M.N_STEPS))
    P4_traces = np.zeros((N_SHOW, M.N_STEPS))
    for i in range(N_SHOW):
        P4_traces[i], P7_traces[i] = run_one(s4, s7, a4, t4, a7, t7, H0,
                                              seed=42 + i)
        print(f"  real {i+1}/{N_SHOW}: P7(2ps)={P7_traces[i,-1]:.3f}  "
              f"max P7={P7_traces[i].max():.3f}")
    P7_mean_20 = P7_traces.mean(0)

    # 500-realisation reference
    mc = np.load(PAPER_DIR / 'exciton_mc.npz')
    P7_full = mc['full_P7']
    # noiseless reference
    P7_noiseless = mc['P7_noiseless']

    # save traces for inspection / later use
    np.savez_compressed(PAPER_DIR / 'mc_realizations.npz',
                        t=M.TLIST, P7=P7_traces, P4=P4_traces,
                        P7_mean_20=P7_mean_20)

    t_fs = M.TLIST * 1e3
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    ax.plot(t_fs, P7_traces.T, color='steelblue', alpha=0.35, lw=0.6)
    ax.plot([], [], color='steelblue', alpha=0.6, lw=1.5,
            label=f'{N_SHOW} individual realisations')
    ax.plot(t_fs, P7_mean_20, color='black', lw=2.2,
            label=f'{N_SHOW}-realisation average')
    ax.plot(t_fs, P7_full, color='#d62728', lw=1.6, ls='--',
            label='500-realisation average')
    ax.plot(t_fs, P7_noiseless, color='gray', lw=1.0, ls=':', alpha=0.7,
            label='noiseless (coherent)')
    ax.set_xlabel('time (fs)')
    ax.set_ylabel(r'population $P_7(t)$')
    ax.set_ylim(-0.05, 1.0); ax.set_xlim(0, t_fs[-1])
    ax.legend(loc='center right')
    # quantify residual oscillation in individual traces
    # (mean number of zero-crossings of d^2P7/dt^2 -- a coherence proxy)
    fig.tight_layout()
    save(fig, 'fig3_realizations.png')

    # ── summary stats ──
    # does each individual trajectory cross P7=0.5 more than once? (coherent)
    crossings = np.array([np.sum(np.abs(np.diff(np.sign(p - 0.5))) > 0)
                          for p in P7_traces])
    print(f"\nindividual P7 traces:")
    print(f"  mean P7(2ps) = {P7_traces[:,-1].mean():.3f} "
          f"(std {P7_traces[:,-1].std():.3f})")
    print(f"  mean max P7  = {P7_traces.max(axis=1).mean():.3f}")
    print(f"  crossings of P7=0.5 per trace: mean={crossings.mean():.1f}  "
          f"range {crossings.min()}-{crossings.max()}")


if __name__ == '__main__':
    main()
