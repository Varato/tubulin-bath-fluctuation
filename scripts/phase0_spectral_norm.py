#!/usr/bin/env python
"""Phase 0 — spectral normalization validation gate.

Goal: confirm PSD computation pipeline + data satisfy ∫ S(f) df == σ².

Two complementary estimators are reported:
  * periodogram  — full-signal FFT with boxcar window. Mathematically exact
                   by Parseval (one-sided). Validates code + data integrity.
  * welch        — Hann-windowed, averaged. This is what Phase 2 will use.
                   Underestimates by a few % due to per-segment demeaning +
                   windowing; the size of that bias is documented here so
                   Phase 2 knows what to expect.

README reports ∫S df / σ² = 0.976-0.999 for the fast dataset (Welch, their
nperseg). PASS criterion here is applied to the periodogram estimate.

Run:  python scripts/phase0_spectral_norm.py
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.signal import welch, periodogram

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load, V_TO_CM, PS_INV_TO_CM, TRP_PRETTY,
                   setup_style, save_fig, trapz, phase_dir)

PASS_MIN, PASS_MAX = 0.998, 1.002   # periodogram: only numerical error

rows = []
peri_ratios = {}
welch_ratios = {}
overall_ok = True
for tag in ['slow', 'fast']:
    d = load(tag)
    dt = d['dt_ps']
    fs = 1.0 / dt                              # ps^-1
    x = d['delta_s_total'] * V_TO_CM           # (N,8) cm^-1
    N = d['n_frames']
    welch_nperseg = min(N, 16384 if tag == 'fast' else 1024)

    peri = np.empty(8)
    welc = np.empty(8)
    print(f"[{tag}] N={N}  dt={dt:.4f} ps  fs={fs:.4g} ps^-1"
          f"  (Nyquist {0.5*fs*PS_INV_TO_CM:.2f} cm^-1)")
    for s in range(8):
        # periodogram (exact): one-sided PSD, demeaned
        f_p, P_p = periodogram(x[:, s], fs=fs, detrend='constant',
                               window='boxcar')
        int_p = float(trapz(P_p, f_p))
        # welch (Phase 2 method)
        f_w, P_w = welch(x[:, s], fs=fs, nperseg=welch_nperseg,
                         detrend='constant')
        int_w = float(trapz(P_w, f_w))

        var = float(x[:, s].var())
        peri[s] = int_p / var
        welc[s] = int_w / var
        rows.append((tag, TRP_PRETTY[s], s, var ** 0.5, peri[s], welc[s]))

    peri_ratios[tag] = peri
    welch_ratios[tag] = welc
    sig = np.sqrt(x.var(0))
    pmin, pmax, pmean = peri.min(), peri.max(), peri.mean()
    wmin, wmax, wmean = welc.min(), welc.max(), welc.mean()
    ok = (pmin >= PASS_MIN) and (pmax <= PASS_MAX)
    overall_ok &= ok
    print(f"   periodogram  ∫Sdf/σ² : min={pmin:.5f}  mean={pmean:.5f}  max={pmax:.5f}"
          f"   ({'PASS' if ok else 'FAIL'} in [{PASS_MIN},{PASS_MAX}])")
    print(f"   welch nperseg={welch_nperseg:<6d} ∫Sdf/σ² : min={wmin:.4f}  mean={wmean:.4f}  max={wmax:.4f}"
          f"   (expected <1: windowing+demean loss)")
    print(f"   σ (cm^-1): min={sig.min():.2f}  mean={sig.mean():.2f}  max={sig.max():.2f}")
    print()

# write table
out = phase_dir(0) / 'spectral_norm.csv'
with out.open('w') as fh:
    fh.write("dataset,site,idx,sigma_cm,periodogram_ratio,welch_ratio\n")
    for r in rows:
        fh.write(f"{r[0]},{r[1]},{r[2]},{r[3]:.4f},{r[4]:.6f},{r[5]:.5f}\n")
print(f"wrote {out}  ({len(rows)} rows)")

# plot
setup_style()
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7, 3.4))
xs = np.arange(8)
w = 0.4
ax.bar(xs - w/2, peri_ratios['fast'], w*0.9, label='fast · periodogram', color='#2c7fb8')
ax.bar(xs + w/2, welch_ratios['fast'], w*0.9, label='fast · welch', color='#41b6c4')
ax.axhline(1.0, color='k', lw=0.8, ls='--', alpha=0.7)
ax.set_xticks(xs); ax.set_xticklabels(TRP_PRETTY, rotation=30, ha='right')
ax.set_ylim(0.6, 1.05)
ax.set_ylabel(r'$\int S(f)\,df\ /\ \sigma^2$')
ax.set_title('Phase 0 — spectral normalization (fast dataset)')
ax.legend(loc='lower right', ncol=2)
fig_path = save_fig(fig, phase_dir(0) / 'spectral_norm.png')
print(f"wrote {fig_path}")

print(f"\nRESULT: {'ALL PASS (periodogram)' if overall_ok else 'FAILURE — STOP'}")
sys.exit(0 if overall_ok else 1)
