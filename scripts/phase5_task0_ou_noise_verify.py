#!/usr/bin/env python
"""Verify: does the SUM of three OU processes give the correct tri-exp ACF?

Single long trajectory (dt=10fs, T=30ns) — resolves all three timescales.
No stitching, no demeaning tricks.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import setup_style, save_fig, phase_dir

AMPS  = [0.50, 0.33, 0.16]
TAUS  = [0.044, 1.70, 2663.0]    # ps
DT    = 0.01                      # ps (10 fs — resolves tau1=44fs)
T_MAX = 30_000.0                  # ps (30 ns — 11x tau3=2.7ns)
N     = int(T_MAX / DT) + 1       # 3,000,001
N_REAL = 5


def gen_ou(n, dt, var, tau, rng):
    decay = np.exp(-dt / tau)
    noise = np.sqrt(var * (1 - decay**2)) * rng.standard_normal(n)
    x = np.empty(n)
    x[0] = np.sqrt(var) * rng.standard_normal()
    for i in range(1, n):
        x[i] = x[i-1] * decay + noise[i]
    return x


def gen_tricolor_sum(n, dt, sigma, amps, taus, rng):
    x = np.zeros(n)
    for A, tau in zip(amps, taus):
        x += gen_ou(n, dt, A * sigma**2, tau, rng)
    return x


def acf_fft(x):
    N = len(x)
    nfft = 1 << int(np.ceil(np.log2(2 * N)))
    F = np.fft.rfft(x - x.mean(), n=nfft)
    ac = np.fft.irfft(F * np.conj(F), n=nfft)[:N]
    return ac / ac[0]


def main():
    P = phase_dir(5)
    setup_style()

    sigma = 998.5   # cm^-1 (Trp4)
    cutoff = N // 4

    print(f"Grid: N={N:,}, dt={DT*1e3:.0f} fs, T={T_MAX/1e3:.0f} ns")
    print(f"ACF reliable to {cutoff*DT/1e3:.1f} ns  (tau3={TAUS[2]/1e3:.1f} ns)")
    print(f"Generating {N_REAL} realizations...")

    rng = np.random.default_rng(42)
    acf_acc = np.zeros(cutoff)
    for i in range(N_REAL):
        x = gen_tricolor_sum(N, DT, sigma, AMPS, TAUS, rng)
        acf_acc += acf_fft(x)[:cutoff]
        print(f"  {i+1}/{N_REAL}  std={x.std():.1f} cm^-1")
    acf_emp = acf_acc / N_REAL

    # Theory
    t = np.arange(cutoff) * DT
    acf_theory = sum(A * np.exp(-t / tau) for A, tau in zip(AMPS, TAUS))

    resid = acf_emp - acf_theory
    rms = np.sqrt(np.mean(resid**2))

    # ── One figure: ACF + residual ──
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9, 7), sharex=True,
        gridspec_kw={'height_ratios': [3, 1]})

    ax1.semilogx(t * 1e3, acf_emp, '-', lw=1.5, color='#333333',
                 alpha=0.8, label='sampled (sum of 3 OU, 5 traj avg)')
    ax1.semilogx(t * 1e3, acf_theory, 'r--', lw=2, alpha=0.7,
                 label=r'theory: $\sum A_k\,e^{-t/\tau_k}$')

    for tau, lbl, c in zip(TAUS,
                            [r'$\tau_1$=44fs', r'$\tau_2$=1.7ps', r'$\tau_3$=2.7ns'],
                            ['#e41a1c', '#377eb8', '#4daf4a']):
        ax1.axvline(tau * 1e3, color=c, ls=':', lw=0.8, alpha=0.5)
        ax1.text(tau * 1e3 * 1.15, 0.45, lbl, fontsize=7, color=c, rotation=90)

    ax1.set_ylabel('C(t) / C(0)')
    ax1.set_title(f'Sum of 3 OU processes vs theory\n'
                  f'dt={DT*1e3:.0f} fs, T={T_MAX/1e3:.0f} ns, '
                  f'{N_REAL} realizations, RMS={rms:.4f}')
    ax1.legend(fontsize=9)
    ax1.set_ylim(-0.05, 1.05)

    ax2.semilogx(t * 1e3, resid, '-', lw=0.6, color='#333333')
    ax2.axhline(0, color='r', lw=0.8, alpha=0.5)
    ax2.set_xlabel('t (fs)')
    ax2.set_ylabel('emp $-$ theory')
    ax2.set_xlim(t[1] * 1e3, t[-1] * 1e3)

    plt.tight_layout()
    save_fig(fig, P / 'ou_sum_acf.png')
    plt.show()

    print(f"\nRMS residual = {rms:.4f}")
    print(f"max |residual| = {np.max(np.abs(resid)):.4f}")
    print(f"Figure: {P}/ou_sum_acf.png")


if __name__ == '__main__':
    main()
