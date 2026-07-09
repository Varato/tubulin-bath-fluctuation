#!/usr/bin/env python
"""Fig A3: verification that the sum of 3 OU processes reproduces the tri-exp ACF.

Loads ou_verify.npz (from run_ou_verify.py). Top: sampled ACF vs theory on log-t.
Bottom: residual. Annotates RMS.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import setup, save, PAPER_DIR


def main():
    setup()
    d = np.load(PAPER_DIR / 'ou_verify.npz')
    t = d['t']; emp = d['acf_emp']; th = d['acf_theory']; rms = float(d['rms'])
    taus = d['taus']
    resid = emp - th

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(11, 6.0), sharex=True,
        gridspec_kw={'height_ratios': [3, 1]})
    ax1.semilogx(t, emp, '-', lw=1.5, color='#333', alpha=0.85,
                 label='sampled (3 OU sum)')
    ax1.semilogx(t, th, 'r--', lw=2, alpha=0.7,
                 label=r'theory $\sum_k f_k e^{-t/T_k}$')
    for tau in taus:
        ax1.axvline(tau, color='gray', ls=':', lw=0.7, alpha=0.4)
    ax1.set_ylabel(r'$C(t)\,/\,C(0)$')
    ax1.set_ylim(-0.05, 1.05)
    ax1.legend(loc='upper right')
    ax1.set_title(f'Sum of three OU processes vs analytic ACF  (RMS = {rms:.3f})')

    ax2.semilogx(t, resid, '-', lw=0.6, color='#333')
    ax2.set_xlabel('lag time $t$ (ps)')
    ax2.set_ylabel('error')
    ax2.set_xlim(t[1], t[-1])
    fig.tight_layout()
    save(fig, 'figA3_ou_verify.png')


if __name__ == '__main__':
    main()
