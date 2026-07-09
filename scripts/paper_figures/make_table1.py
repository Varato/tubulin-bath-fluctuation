#!/usr/bin/env python
"""Table 1: per-site fluctuation characterization.

site | sigma_total | sigma/J | f1 | T1 | f2 | T2 | f3 | T3   (+ mean row)
Data: sigma_matrix.csv (slow) + acf_fit_triexp.csv. Writes table1.md + .csv.
"""
from __future__ import annotations
import sys, csv
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import TRP_NAMES, PAPER_DIR
from utils import J_COUPLING

P1 = PAPER_DIR.parent / 'phase1_basic_stats'
P2 = PAPER_DIR.parent / 'phase2_timescale'


def main():
    sigma = {}
    with open(P1 / 'sigma_matrix.csv') as f:
        for r in csv.DictReader(f):
            if r['dataset'] == 'slow':
                sigma[int(r['idx'])] = float(r['sigma_total_cm'])
    tri = {}
    R2 = {}
    # corrected tri-exp: tau3 anchored to the unbiased slow-only estimate
    # (refit_triexp_corrected.py) — single source for table, figure, and R^2.
    with open(PAPER_DIR / 'corrected_triexp.csv') as f:
        for i, r in enumerate(csv.DictReader(f)):
            tri[i] = [float(r['A1']), float(r['tau1_ps']),
                      float(r['A2']), float(r['tau2_ps']),
                      float(r['A3']), float(r['tau3_ps'])]
            R2[i] = float(r['R2'])

    rows = []
    for i in range(8):
        s = sigma[i]; a1, t1, a2, t2, a3, t3 = tri[i]
        rows.append([TRP_NAMES[i], s, s / J_COUPLING,
                     a1, t1, a2, t2, a3, t3])
    # mean row
    arr = np.array([[r[j] for j in range(1, 9)] for r in rows])
    mean = arr.mean(axis=0)
    rows.append(['mean', *mean])

    # sanity: amplitudes sum to 1
    for i in range(8):
        assert abs(tri[i][0] + tri[i][2] + tri[i][4] - 1.0) < 2e-3  # CSV rounding
    print("amplitudes sum to 1 per site (within CSV rounding): OK")

    # csv
    hdr = 'site,sigma_total_cm,sigma_over_J,f1,T1_ps,f2,T2_ps,f3,T3_ps'
    with open(PAPER_DIR / 'table1.csv', 'w') as f:
        f.write(hdr + '\n')
        for r in rows:
            f.write(','.join(f'{v:.4g}' if isinstance(v, float) else str(v)
                             for v in r) + '\n')

    # markdown — unified time unit: ps
    def fmt_tau(t):
        if t < 0.1:
            return f'{t:.3f}'      # T1: ~0.03–0.07 ps
        if t < 10:
            return f'{t:.2f}'      # T2: ~0.4–9 ps
        return f'{t:.0f}'          # T3: ~450–8500 ps
    lines = ['| site | σ (cm⁻¹) | σ/J | f₁ | T₁ (ps) | f₂ | T₂ (ps) | f₃ | T₃ (ps) |',
             '|---|---|---|---|---|---|---|---|---|']
    for r in rows:
        nm = r[0]
        if nm == 'mean':
            lines.append(f'| **mean** | **{r[1]:.0f}** | **{r[2]:.1f}** | '
                         f'{r[3]:.2f} | {fmt_tau(r[4])} | {r[5]:.2f} | '
                         f'{fmt_tau(r[6])} | {r[7]:.2f} | {fmt_tau(r[8])} |')
        else:
            lines.append(f'| {nm} | {r[1]:.0f} | {r[2]:.1f} | '
                         f'{r[3]:.2f} | {fmt_tau(r[4])} | {r[5]:.2f} | '
                         f'{fmt_tau(r[6])} | {r[7]:.2f} | {fmt_tau(r[8])} |')
    md = '\n'.join(lines) + '\n'
    (PAPER_DIR / 'table1.md').write_text(md)
    print(f"saved {PAPER_DIR/'table1.md'}")
    print(f"saved {PAPER_DIR/'table1.csv'}")
    print(md)


if __name__ == '__main__':
    main()
