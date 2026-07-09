#!/usr/bin/env python
"""Compute OU-sum verification: does sum of 3 OU processes reproduce the
tri-exp ACF? Single long trajectory (dt=10fs, T=30ns). Saves ou_verify.npz.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import PAPER_DIR

AMPS = [0.51, 0.33, 0.16]
TAUS = [0.046, 1.77, 2663.0]      # ps (corrected system means)
DT = 0.01                          # ps
T_MAX = 30_000.0                   # ps (30 ns)
N = int(T_MAX / DT) + 1
N_REAL = 5


def gen_ou(n, dt, var, tau, rng):
    decay = np.exp(-dt / tau)
    noise = np.sqrt(var * (1 - decay**2)) * rng.standard_normal(n)
    x = np.empty(n)
    x[0] = np.sqrt(var) * rng.standard_normal()
    for i in range(1, n):
        x[i] = x[i - 1] * decay + noise[i]
    return x


def acf_fft(x):
    N = len(x)
    nfft = 1 << int(np.ceil(np.log2(2 * N)))
    F = np.fft.rfft(x - x.mean(), n=nfft)
    ac = np.fft.irfft(F * np.conj(F), n=nfft)[:N]
    return ac / ac[0]


def main():
    sigma = 998.5   # cm^-1 (Trp4), arbitrary — ACF is normalized
    cutoff = N // 4
    print(f"N={N:,} dt={DT*1e3:.0f}fs T={T_MAX/1e3:.0f}ns  cutoff={cutoff*DT/1e3:.1f}ns")
    rng = np.random.default_rng(42)
    acc = np.zeros(cutoff)
    for i in range(N_REAL):
        x = np.zeros(N)
        for A, tau in zip(AMPS, TAUS):
            x += gen_ou(N, DT, A * sigma**2, tau, rng)
        acc += acf_fft(x)[:cutoff]
        print(f"  {i+1}/{N_REAL} std={x.std():.1f}")
    acf_emp = acc / N_REAL
    t = np.arange(cutoff) * DT
    acf_theory = sum(A * np.exp(-t / tau) for A, tau in zip(AMPS, TAUS))
    rms = float(np.sqrt(np.mean((acf_emp - acf_theory)**2)))
    np.savez_compressed(PAPER_DIR / 'ou_verify.npz',
                        t=t, acf_emp=acf_emp, acf_theory=acf_theory,
                        rms=rms, amps=np.array(AMPS), taus=np.array(TAUS))
    print(f"RMS residual = {rms:.4f}")
    print(f"saved {PAPER_DIR/'ou_verify.npz'}")


if __name__ == '__main__':
    main()
