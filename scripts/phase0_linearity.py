#!/usr/bin/env python
"""Phase 0 — linearity validation gate.

Verifies E_total == E_protein + E_water + E_nucleotide + E_ions to numerical
precision on BOTH datasets. RESEARCH_PLAN §1.6 lists this as a prerequisite;
README reports max ratio ~ 4e-14.

If this check fails the source-decomposition in every later phase is invalid.

Run:  python scripts/phase0_linearity.py
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import load, FILES, TRP_PRETTY, COMPONENTS, phase_dir

PASS_THRESHOLD = 1e-9   # numerical precision; README saw 4e-14

rows = []
overall_ok = True
for tag in ['slow', 'fast']:
    d = load(tag)
    E_sum = sum(d[f'E_{c}'] for c in COMPONENTS)         # (N,8,3)
    resid = d['E_total'] - E_sum
    denom = np.abs(d['E_total']).max()
    if denom == 0:
        print(f"[{tag}] FATAL: max|E_total| == 0"); sys.exit(2)
    max_ratio = float(np.abs(resid).max() / denom)
    rms_ratio = float(np.sqrt((resid ** 2).mean()) / denom)

    ok = max_ratio < PASS_THRESHOLD
    overall_ok &= ok
    print(f"[{tag}] N={d['n_frames']:>7d}  dt={d['dt_ps']:>7.3f} ps")
    print(f"   max|E_total - Σ E_sources| / max|E_total|  = {max_ratio:.3e}"
          f"   ({'PASS' if ok else 'FAIL'} < {PASS_THRESHOLD:.0e})")
    print(f"   rms|E_total - Σ E_sources| / max|E_total|  = {rms_ratio:.3e}")

    # per-site max ratio for the CSV
    per_site = (np.abs(resid).reshape(d['n_frames'], 8, 3).max(axis=(0, 2))
                / np.abs(d['E_total']).reshape(d['n_frames'], 8, 3).max(axis=(0, 2)))
    for s, r in enumerate(per_site):
        rows.append((tag, TRP_PRETTY[s], s, float(r)))
    print()

# write per-site table
out = phase_dir(0) / 'linearity.csv'
with out.open('w') as fh:
    fh.write("dataset,site,idx,max_ratio\n")
    for r in rows:
        fh.write(f"{r[0]},{r[1]},{r[2]},{r[3]:.6e}\n")
print(f"wrote {out}  ({len(rows)} rows)")

print(f"\nRESULT: {'ALL PASS' if overall_ok else 'FAILURE — STOP'}")
sys.exit(0 if overall_ok else 1)
