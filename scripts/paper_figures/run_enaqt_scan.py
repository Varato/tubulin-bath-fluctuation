#!/usr/bin/env python
"""ENAQT parameter scan: T3 (fixed at full MD) + additional OU (variable).

Two orthogonal scans decomposing the ENAQT parameter space:
  (b) Fix σ_OU, scan T_OU  — ENAQT vs bath correlation time.
  (d) Fix T_OU, scan σ_OU  — ENAQT vs bath amplitude.

T3 noise is always at full MD strength (Anderson localization baseline).
The additional OU is a generic controllable bath on top.

Expected curves (both rise-then-fall):
  T-scan:  T→∞ (static) → localization; T≈1/J (resonance crossings) → peak;
           T→0 (white-noise limit, σ²T→0) → no effect → back to T3 baseline.
  σ-scan: σ→0 → only T3; σ≈J → ENAQT peak; σ≫J → overdamped suppression.

Output: results/paper_figures/enaqt_scan.npz
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import CM_TO_RADPS, PAPER_DIR
import run_exciton_mc as M

N_REAL = 500


def gen_t3_plus_ou(n_steps, dt, sigma_m, f3_m, T3_m, sigma_ou, T_ou, rng):
    """T3 component + additional OU for one site. Returns noise array (cm^-1)."""
    noise = M.gen_ou(n_steps, dt, f3_m * sigma_m**2, T3_m, rng)
    if sigma_ou > 0 and T_ou > 0:
        noise += M.gen_ou(n_steps, dt, sigma_ou**2, T_ou, rng)
    return noise


def run_mc(sig, tri, sigma_ou, T_ou, H0, n_real):
    """MC with T3 + additional OU. Returns P7(t=4ps) mean and SEM."""
    s4, s7 = sig[M.I4], sig[M.I7]
    a4, t4 = tri[M.I4]; a7, t7 = tri[M.I7]
    f3_4, T3_4 = a4[2], t4[2]
    f3_7, T3_7 = a7[2], t7[2]
    psi0 = np.array([1., 0.], dtype=complex)
    P7_end = np.zeros(n_real)
    for i in range(n_real):
        rng = np.random.default_rng(42 + i)
        n4 = gen_t3_plus_ou(M.N_STEPS, M.DT, s4, f3_4, T3_4,
                            sigma_ou, T_ou, rng)
        n7 = gen_t3_plus_ou(M.N_STEPS, M.DT, s7, f3_7, T3_7,
                            sigma_ou, T_ou, rng)
        psi = psi0.copy()
        for j in range(1, M.N_STEPS):
            H = H0 + np.diag([n4[j], n7[j]]) * CM_TO_RADPS
            psi = expm(-1j * H * M.DT) @ psi
        P7_end[i] = abs(psi[1])**2
    return P7_end.mean(), P7_end.std() / np.sqrt(n_real)


def main():
    sig, tri = M.load_site_params()
    J  = M.J_CM   * CM_TO_RADPS
    dE = M.D_E_CM * CM_TO_RADPS
    H0 = J * np.array([[0, 1], [1, 0]], dtype=complex) \
       + (dE / 2) * np.array([[1, 0], [0, -1]], dtype=complex)

    print(f"T3+OU ENAQT scan. J={M.J_CM} cm^-1, dE={M.D_E_CM} cm^-1")
    print(f"Grid: {M.N_STEPS} steps, dt={M.DT*1e3:.0f}fs, T={M.T_MAX}ps, N_real={N_REAL}\n")

    # ── (b) T-scan: fix σ_OU, scan T_OU ──
    SIGMA_FIX = 300.0   # cm^-1, the fixed additional-OU amplitude
    T_values = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05,
                0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0]
    print(f"=== (b) T-scan (σ_OU={SIGMA_FIX} cm^-1 fixed) ===")
    t_scan_p7 = []; t_scan_sem = []
    t0 = time.time()
    for T_ou in T_values:
        p7, sem = run_mc(sig, tri, SIGMA_FIX, T_ou, H0, N_REAL)
        t_scan_p7.append(p7); t_scan_sem.append(sem)
        print(f"  T={T_ou:8.3f} ps  P7={p7:.3f} +- {sem:.3f}  "
              f"({(time.time()-t0):.0f}s)")
    t_scan_p7 = np.array(t_scan_p7); t_scan_sem = np.array(t_scan_sem)

    # ── (d) σ-scan: fix T_OU, scan σ_OU ──
    T_FIX = 0.05   # ps, the fixed additional-OU correlation time (≈T1≈1/J)
    SIGMA_values = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500,
                    1000, 2000, 5000]
    print(f"\n=== (d) σ-scan (T_OU={T_FIX} ps fixed) ===")
    s_scan_p7 = []; s_scan_sem = []
    t0 = time.time()
    for sigma_ou in SIGMA_values:
        p7, sem = run_mc(sig, tri, sigma_ou, T_FIX, H0, N_REAL)
        s_scan_p7.append(p7); s_scan_sem.append(sem)
        print(f"  σ={sigma_ou:6.0f} cm^-1  P7={p7:.3f} +- {sem:.3f}  "
              f"({(time.time()-t0):.0f}s)")
    s_scan_p7 = np.array(s_scan_p7); s_scan_sem = np.array(s_scan_sem)

    np.savez_compressed(
        PAPER_DIR / 'enaqt_scan.npz',
        # T-scan
        sigma_fix=SIGMA_FIX, T_values=np.array(T_values),
        T_scan_P7=t_scan_p7, T_scan_sem=t_scan_sem,
        # sigma-scan
        T_fix=T_FIX, sigma_values=np.array(SIGMA_values),
        sigma_scan_P7=s_scan_p7, sigma_scan_sem=s_scan_sem,
        # reference
        J_cm=M.J_CM,
    )
    print(f"\nsaved -> {PAPER_DIR / 'enaqt_scan.npz'}")
    print(f"\nSummary:")
    print(f"  T-scan peak: P7={t_scan_p7.max():.3f} at T={T_values[t_scan_p7.argmax()]:.3f} ps")
    print(f"  σ-scan peak: P7={s_scan_p7.max():.3f} at σ={SIGMA_values[s_scan_p7.argmax()]:.0f} cm^-1")


if __name__ == '__main__':
    main()
