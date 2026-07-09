#!/usr/bin/env python
"""Local angular deviation of each Trp transition dipole.

Identical methodology to:
  microtubule-quantum-energy-transfer/13_indole_mobility/compute_local_ang_dev.py

Adapted to read from our NPZ pipeline (utils.load) instead of the JSON
fluctuations file.  The slow trajectory (10-50 ns, 4001 frames, dt=10 ps)
matches the reference analysis window.

Method:
  1. Extract lab-frame dipole directions n̂_m(t) from dmu array (normalised).
  2. Per-frame Procrustes alignment removes global protein rotation R(t):
       solve Wahba's problem  min_R Σ_m |R⟨n̂_m⟩ − n̂_m(t)|²  via SVD.
  3. Compute mean angular tilt θ_loc = ⟨arccos(n̂_corr(t)·⟨n̂_corr⟩)⟩.

θ_loc is the orientational analogue of RMSF — how many degrees the dipole
typically tilts from its average orientation, after removing global rotation.
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (TRP_PRETTY, TRP_LABELS, load, setup_style,
                   write_csv, RESULTS)


def procrustes_rotation(ref_dirs, cur_dirs):
    """Optimal rotation R mapping ref_dirs → cur_dirs (Wahba's problem).

    Given unit vectors a_m = ref_dirs[m] and b_m = cur_dirs[m],
    finds R that minimises  Σ_m |R a_m − b_m|².

    Returns R (3×3) such that  cur_dirs ≈ (R @ ref_dirs.T).T.
    """
    H = ref_dirs.T @ cur_dirs        # (3,3) cross-covariance
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    return Vt.T @ D @ U.T


def angular_deviation(directions):
    """Mean angular tilt of a set of vectors from their mean direction.

    directions : (N, 3) array.  Need NOT be unit vectors — they are
    normalised internally before computing the angle.
    Returns mean angle in degrees.
    """
    normed = directions / np.linalg.norm(directions, axis=1, keepdims=True)
    mean_dir = normed.mean(axis=0)
    mean_dir /= np.linalg.norm(mean_dir)
    cos_ang = np.clip(normed @ mean_dir, -1, 1)
    return np.degrees(np.arccos(cos_ang)).mean()


def main():
    setup_style()
    out = RESULTS / 'structural_correlation'
    out.mkdir(parents=True, exist_ok=True)

    # --- Load slow trajectory ---
    d = load('slow')
    N = d['n_frames']
    dt = d['dt_ps']
    print(f"Loaded slow trajectory: {N} frames, dt={dt} ps, "
          f"T={N*dt/1000:.1f} ns")

    # --- Extract lab-frame dipole directions: shape (8, N, 3) ---
    dmu = d['dmu']                  # (N, 8, 3) raw dipole vectors
    dmu_lab = np.transpose(dmu, (1, 0, 2))  # (8, N, 3)
    # Normalise to unit vectors (direction only; magnitude carries oscillator
    # strength which we don't want weighting the Procrustes fit)
    dmu_lab /= np.linalg.norm(dmu_lab, axis=2, keepdims=True)

    # --- Reference: time-averaged directions per Trp ---
    ref_dirs = dmu_lab.mean(axis=1)   # (8, 3)

    # --- Per-frame Procrustes correction ---
    dmu_corrected = np.zeros_like(dmu_lab)
    for k in range(N):
        cur = dmu_lab[:, k, :]           # (8, 3)
        R = procrustes_rotation(ref_dirs, cur)
        dmu_corrected[:, k, :] = (R.T @ cur.T).T
    # Re-normalise
    dmu_corrected /= np.linalg.norm(dmu_corrected, axis=2, keepdims=True)

    # --- Compute local angular deviation ---
    lab_ang = np.array([angular_deviation(dmu_lab[i])       for i in range(8)])
    cor_ang = np.array([angular_deviation(dmu_corrected[i]) for i in range(8)])

    # --- Diagnostic: mean-direction deviation (should drop after correction) ---
    lab_mean_dev = angular_deviation(dmu_lab.mean(axis=0))
    cor_mean_dev = angular_deviation(dmu_corrected.mean(axis=0))

    # --- Print ---
    print(f"\n{'site':<8} {'lab (deg)':>10} {'local (deg)':>12}")
    print("-" * 35)
    for i in range(8):
        print(f"{TRP_PRETTY[i]:<8} {lab_ang[i]:>10.1f} {cor_ang[i]:>12.1f}")
    print(f"\nGlobal rotation check (mean-direction deviation):")
    print(f"  Lab frame:        {lab_mean_dev:.1f}°")
    print(f"  After Procrustes: {cor_mean_dev:.1f}°  (residual)")

    # --- Save JSON ---
    out_json = {
        "description": "Local angular deviation of Trp transition dipole, "
                       "Procrustes-corrected for global protein rotation.",
        "source": "slow trajectory (10-50 ns, 4001 frames, dt=10 ps)",
        "n_frames": int(N),
        "trp_labels": TRP_LABELS,
        "trp_pretty": TRP_PRETTY,
        "angular_deviation_lab_deg": lab_ang.tolist(),
        "angular_deviation_local_deg": cor_ang.tolist(),
        "global_rotation_residual_deg": float(cor_mean_dev),
    }
    json_path = out / 'indole_mobility.json'
    with open(json_path, 'w') as f:
        json.dump(out_json, f, indent=2)
    print(f"\nSaved: {json_path}")

    # --- Save CSV ---
    rows = [(TRP_PRETTY[i], f'{lab_ang[i]:.2f}', f'{cor_ang[i]:.2f}')
            for i in range(8)]
    write_csv(out / 'indole_mobility.csv',
              'site,ang_dev_lab_deg,ang_dev_local_deg', rows)
    print(f"Saved: {out}/indole_mobility.csv")


if __name__ == '__main__':
    main()
