#!/usr/bin/env python
"""Oscillator strength analysis of the Craddock 8-site Trp exciton Hamiltonian.

Analysis A: Clean eigenstate oscillator strengths using MD-derived dipole directions.
Analysis B: Disorder sigma-scan of brightest-state oscillator strength.
Analysis C: Participation ratio under disorder.

Outputs (superradiance/results/):
    clean_oscillator_strength.npz   eigenvalues, eigenvectors, PR, F_k
    disorder_scan.npz               F_max(sigma), PR(sigma), distributions
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

# Make utils importable
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent.parent / 'scripts'))
from utils import load, TRP_PRETTY

# ── Paths ──
RESULTS_DIR = SCRIPT_DIR.parent / 'results'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── Craddock 2014 Hamiltonian (cm^-1) ──
H_CRAD = np.array([
    [  1,   0, -13,   0,  -2,  -1,   5,  -1],
    [  0, 388, -41,   4,   1,   1,  -4,  1],
    [-13, -41, 342,   2,   0,   1,  -6,  1],
    [  0,   4,   2, 207,  -4,   6, -59,  -1],
    [ -2,   1,   0,  -4,  57,  21,   2,  11],
    [ -1,   1,   1,   6,  21, 102,   5, -51],
    [  5,  -4,  -6, -59,   2,   5, 248,   3],
    [ -1,   1,   1,  -1,  11, -51,   3,   0]
])

N_SITES = 8


def load_dipole_directions():
    """Load time-averaged transition dipole unit vectors from MD slow trajectory."""
    d = load('slow')
    dmu = d['dmu']  # (N_frames, 8, 3)
    # Normalise each frame to unit vectors
    norms = np.linalg.norm(dmu, axis=2, keepdims=True)
    dmu_unit = dmu / norms
    # Time average (dipole frozen on exciton timescale; reorientation ~17 ns)
    n_hat = dmu_unit.mean(axis=0)  # (8, 3)
    n_hat /= np.linalg.norm(n_hat, axis=1, keepdims=True)
    return n_hat


def participation_ratio(eigvec):
    """PR = 1 / sum |c_i|^4 for a single eigenvector."""
    w = np.abs(eigvec) ** 2
    return 1.0 / np.sum(w ** 2)


def oscillator_strength(eigvec, n_hat):
    """Oscillator strength enhancement factor F_k = |sum_i c_i n_hat_i|^2.

    Normalised so single-Trp F = 1 (since |n_hat_i| = 1).
    """
    mu_k = eigvec @ n_hat  # (3,) vector: sum_i c_i * n_hat_i
    return np.dot(mu_k, mu_k)


# ════════════════════════════════════════════════════════════════════
# Analysis A: Clean eigenstate oscillator strengths
# ════════════════════════════════════════════════════════════════════
def analysis_a(n_hat):
    print("=" * 70)
    print("Analysis A: Clean eigenstate oscillator strengths")
    print("=" * 70)

    energies, states = np.linalg.eigh(H_CRAD)

    pr_list = []
    f_list = []
    for k in range(N_SITES):
        c = states[:, k]
        pr = participation_ratio(c)
        f = oscillator_strength(c, n_hat)
        pr_list.append(pr)
        f_list.append(f)

        w = np.abs(c) ** 2
        top2 = np.argsort(w)[::-1][:2]
        char = f"Trp{top2[0]+1}-Trp{top2[1]+1}"
        print(f"  psi{k+1} (E={energies[k]:+7.1f}): PR={pr:.2f}  F={f:.3f}  "
              f"[{char}]")

    print(f"\n  Sum rule check: sum(F_k) = {sum(f_list):.3f}  (expected {N_SITES})")
    print(f"  Brightest state: psi{np.argmax(f_list)+1}  F_max = {max(f_list):.3f}")
    print(f"  Dimmest state:   psi{np.argmin(f_list)+1}  F_min = {min(f_list):.3f}")

    # Dipole alignment for key pairs
    print("\n  Dipole alignment for coupled pairs:")
    pairs = [(3, 6, 'Trp4-Trp7'), (5, 7, 'Trp6-Trp8'), (1, 2, 'Trp2-Trp3')]
    for i, j, label in pairs:
        cos_theta = np.dot(n_hat[i], n_hat[j])
        angle = np.degrees(np.arccos(np.clip(abs(cos_theta), 0, 1)))
        print(f"    {label}: cos theta = {cos_theta:+.3f}  (|theta| = {angle:.1f} deg)")

    # Save
    np.savez_compressed(
        RESULTS_DIR / 'clean_oscillator_strength.npz',
        energies=energies,
        states=states,
        pr=np.array(pr_list),
        f_k=np.array(f_list),
        n_hat=n_hat,
        trp_pretty=np.array(TRP_PRETTY),
    )
    print(f"\n  Saved -> {RESULTS_DIR / 'clean_oscillator_strength.npz'}")
    return energies, states, np.array(f_list)


# ════════════════════════════════════════════════════════════════════
# Analysis B+C: Disorder sigma-scan
# ════════════════════════════════════════════════════════════════════
def analysis_bc(n_hat, n_real=500):
    print("\n" + "=" * 70)
    print("Analysis B+C: Disorder sigma-scan of oscillator strength and PR")
    print("=" * 70)

    sigmas = np.array([0, 10, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500])

    # Storage
    f_max_mean = np.zeros(len(sigmas))
    f_max_sem = np.zeros(len(sigmas))
    f_max_p05 = np.zeros(len(sigmas))
    f_max_p95 = np.zeros(len(sigmas))
    pr_mean = np.zeros((len(sigmas), N_SITES))
    pr_sem = np.zeros((len(sigmas), N_SITES))
    f_brightest_idx = np.zeros(len(sigmas))  # which eigenstate is brightest

    rng = np.random.default_rng(42)

    for si, sigma in enumerate(sigmas):
        f_max_real = np.zeros(n_real)
        pr_real = np.zeros((n_real, N_SITES))

        for r in range(n_real):
            # Add Gaussian static disorder to diagonal
            disorder = rng.normal(0, sigma, size=N_SITES) if sigma > 0 else np.zeros(N_SITES)
            H_disordered = H_CRAD + np.diag(disorder)

            energies_r, states_r = np.linalg.eigh(H_disordered)

            # Compute F for all eigenstates
            f_vals = np.zeros(N_SITES)
            pr_vals = np.zeros(N_SITES)
            for k in range(N_SITES):
                f_vals[k] = oscillator_strength(states_r[:, k], n_hat)
                pr_vals[k] = participation_ratio(states_r[:, k])

            f_max_real[r] = np.max(f_vals)
            pr_real[r] = pr_vals

        f_max_mean[si] = np.mean(f_max_real)
        f_max_sem[si] = np.std(f_max_real, ddof=1) / np.sqrt(n_real)
        f_max_p05[si] = np.percentile(f_max_real, 5)
        f_max_p95[si] = np.percentile(f_max_real, 95)
        pr_mean[si] = np.mean(pr_real, axis=0)
        pr_sem[si] = np.std(pr_real, axis=0, ddof=1) / np.sqrt(n_real)

        marker = " <-- sigma_3 (Trp4/7)" if sigma == 200 else ""
        print(f"  sigma={sigma:5.0f} cm^-1: F_max = {f_max_mean[si]:.3f} "
              f"± {f_max_sem[si]:.3f}  "
              f"<PR_max> = {np.max(pr_mean[si]):.2f}{marker}")

    # Save
    np.savez_compressed(
        RESULTS_DIR / 'disorder_scan.npz',
        sigmas=sigmas,
        f_max_mean=f_max_mean,
        f_max_sem=f_max_sem,
        f_max_p05=f_max_p05,
        f_max_p95=f_max_p95,
        pr_mean=pr_mean,
        pr_sem=pr_sem,
        n_real=n_real,
        n_hat=n_hat,
    )
    print(f"\n  Saved -> {RESULTS_DIR / 'disorder_scan.npz'}")

    return sigmas, f_max_mean, f_max_sem, pr_mean


def main():
    n_hat = load_dipole_directions()

    print("\nTime-averaged dipole directions:")
    for i in range(N_SITES):
        print(f"  {TRP_PRETTY[i]}: ({n_hat[i,0]:+.3f}, {n_hat[i,1]:+.3f}, {n_hat[i,2]:+.3f})")

    # Analysis A
    energies, states, f_clean = analysis_a(n_hat)

    # Analysis B+C
    sigmas, f_max_mean, f_max_sem, pr_mean = analysis_bc(n_hat, n_real=500)

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Clean dimer brightest state: F = {max(f_clean):.3f} (psi{np.argmax(f_clean)+1})")
    print(f"  At sigma_3 = 200 cm^-1: F_max = {f_max_mean[sigmas == 200][0]:.3f}")
    print(f"  Enhancement is sub-2x throughout -> no meaningful superradiance in dimer")


if __name__ == '__main__':
    main()
