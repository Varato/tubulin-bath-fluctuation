#!/usr/bin/env python
"""Refit the per-site tri-exp with the slow timescale anchored to the
unbiased slow-only estimate (acf_fit_separate.csv). This removes the bias
that the 10 ps stitch discontinuity introduces into tau3 in the stitched
tri-exp (acf_fit_triexp.csv).

Model per site:
    C(t) = A1 exp(-t/tau1) + A2 exp(-t/tau2) + A3 exp(-t/tau3)
with tau3 FIXED = tau_slow (slow-only fit), A3 = 1 - A1 - A2, and
A1, tau1, A2, tau2 floated on the stitched ACF.

Saves corrected_triexp.csv (A1,tau1,A2,tau2,A3,tau3,R2 per site).
"""
from __future__ import annotations
import sys, csv
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import PAPER_DIR

P2 = PAPER_DIR.parent / 'phase2_timescale'


def load():
    d = np.load(P2 / 'acf_stitched.npz')
    t = d['t_lag_ps']; acf = d['acf']                      # (N,8)
    # stitched tri-exp (initial guesses)
    st = {}
    for r in csv.DictReader(open(P2 / 'acf_fit_triexp.csv')):
        i = int(r['idx'])
        st[i] = ([float(r['A1']), float(r['A2']), float(r['A3'])],
                 [float(r['tau1_ps']), float(r['tau2_ps']), float(r['tau3_ps'])])
    # slow-only tau_slow (anchor)
    slow = {}
    for r in csv.DictReader(open(P2 / 'acf_fit_separate.csv')):
        slow[int(r['idx'])] = float(r['tau_slow_slowfit_ps'])
    return t, acf, st, slow


def fit_one(t, y, tau3_fixed, p0):
    A1_0, tau1_0, A2_0, tau2_0 = p0
    def model(t, A1, tau1, A2, tau2):
        A3 = 1.0 - A1 - A2
        return (A1 * np.exp(-t / tau1) + A2 * np.exp(-t / tau2)
                + A3 * np.exp(-t / tau3_fixed))
    lo = [0.0, 1e-4, 0.0, 1e-4]
    hi = [1.0, 10.0, 1.0, 100.0]   # tau in ps; tau1,tau2 < 10 ps (sub-stitch)
    popt, _ = curve_fit(model, t, y, p0=[A1_0, tau1_0, A2_0, tau2_0],
                        bounds=(lo, hi), maxfev=20000)
    fit = model(t, *popt)
    ss_res = float(np.sum((y - fit) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot
    A1, tau1, A2, tau2 = popt
    A3 = 1.0 - A1 - A2
    return A1, tau1, A2, tau2, A3, tau3_fixed, r2


def main():
    t, acf, st, slow = load()
    rows = []
    print(f"{'site':<6}{'A1':>7}{'T1':>8}{'A2':>7}{'T2':>8}{'A3':>7}{'T3(ns)':>9}{'R2':>8}")
    print('-' * 60)
    for i in range(8):
        A, tau = st[i]
        p0 = [A[0], tau[0], A[1], tau[1]]
        A1, t1, A2, t2, A3, t3, r2 = fit_one(t, acf[:, i], slow[i], p0)
        name = ['Trp1','Trp2','Trp3','Trp4','Trp5','Trp6','Trp7','Trp8'][i]
        print(f"{name:<6}{A1:>7.3f}{t1*1000:>5.0f}fs{A2:>7.3f}{t2*1000:>5.0f}fs"
              f"{A3:>7.3f}{t3/1000:>9.2f}{r2:>8.4f}")
        rows.append([name, A1, t1, A2, t2, A3, t3, r2])
    # save csv
    with open(PAPER_DIR / 'corrected_triexp.csv', 'w') as f:
        f.write('site,A1,tau1_ps,A2,tau2_ps,A3,tau3_ps,R2\n')
        for r in rows:
            f.write(f'{r[0]},{r[1]:.4f},{r[2]:.5f},{r[3]:.4f},{r[4]:.5f},'
                    f'{r[5]:.4f},{r[6]:.3f},{r[7]:.5f}\n')
    R2 = np.array([r[7] for r in rows])
    t3 = np.array([r[6] for r in rows])
    A3 = np.array([r[5] for r in rows])
    print(f"\nmean R2 = {R2.mean():.4f} (min {R2.min():.4f})")
    print(f"mean T3 = {t3.mean():.0f} ps = {t3.mean()/1000:.2f} ns")
    print(f"mean A3 = {A3.mean():.3f}")
    print(f"saved {PAPER_DIR/'corrected_triexp.csv'}")


if __name__ == '__main__':
    main()
