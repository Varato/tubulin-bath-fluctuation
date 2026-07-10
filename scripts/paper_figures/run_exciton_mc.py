#!/usr/bin/env python
"""Corrected exciton MC for the Trp4-Trp7 (alphaW407-betaW344) dimer.

Fixes the sigma mismatch in phase5_task1_exciton_mc.py (which used betaW397's
sigma). Here we use the physically correct partner of Trp4: Trp7 = betaW344,
the only site with strong coupling to Trp4 (J = -59 cm^-1; all others |J|<=6).

  H(t) = J sx + (dE/2) sz + ds4(t)|0><0| + ds7(t)|1><1|

with per-site sigma and tri-exp OU noise from Phase 1/2 (slow trajectory), and
detuning dE = 41 cm^-1 from the Craddock 2014 diagonal (207 vs 248).

Also computes the Craddock 2014 Haken-Strobl Lindblad baseline (gamma = 50 cm^-1)
on (a) the same 2-site dimer and (b) the full 8-site Hamiltonian.

Outputs (results/paper_figures/):
    exciton_mc.npz         dimer trajectory MC, all configs + noiseless
    lindblad_2site.npz     2-site Haken-Strobl, gamma=50
    lindblad_8site.npz     full 8-site Craddock, gamma=50
"""
from __future__ import annotations
import sys, csv
from pathlib import Path
import numpy as np
import qutip as qt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import CM_TO_RADPS, PAPER_DIR

# ── Craddock 2014 Hamiltonian (cm^-1), diag = site energies ──
H_CRAD = np.array([
    [  1,   0, -13,   0,  -2,  -1,   5,  -1],
    [  0, 388, -41,   4,   1,   1,  -4,   1],
    [-13, -41, 342,   2,   0,   1,  -6,   1],
    [  0,   4,   2, 207,  -4,   6, -59,  -1],
    [ -2,   1,   0,  -4,  57,  21,   2,  11],
    [ -1,   1,   1,   6,  21, 102,   5, -51],
    [  5,  -4,  -6, -59,   2,   5, 248,   3],
    [ -1,   1,   1,  -1,  11, -51,   3,   0]
])

I4, I7 = 3, 6          # Trp4 = alphaW407, Trp7 = betaW344
J_CM   = H_CRAD[I4, I7]                 # -59
D_E_CM = H_CRAD[I4, I4] - H_CRAD[I7, I7]  # 207-248 = -41
GAMMA_CM = 50.0                          # Craddock 2014 dephasing

# ── load per-site sigma (slow traj) + tri-exp from Phase 1/2 CSVs ──
def load_site_params():
    sig = {}
    with open(PAPER_DIR.parent / 'phase1_basic_stats' / 'sigma_matrix.csv') as f:
        for r in csv.DictReader(f):
            if r['dataset'] == 'slow':
                sig[int(r['idx'])] = float(r['sigma_total_cm'])
    tri = {}
    # corrected tri-exp (tau3 anchored to slow-only) — consistent with Table 1
    with open(PAPER_DIR / 'corrected_triexp.csv') as f:
        for i, r in enumerate(csv.DictReader(f)):
            tri[i] = (np.array([float(r['A1']), float(r['A2']), float(r['A3'])]),
                      np.array([float(r['tau1_ps']), float(r['tau2_ps']),
                                float(r['tau3_ps'])]))
    return sig, tri

# ── grid ──
DT, T_MAX = 0.002, 4.0
TLIST = np.arange(0, T_MAX + DT / 2, DT)
N_STEPS = len(TLIST)
N_REAL = 500

# ── 2-site operators ──
k0, k1 = qt.basis(2, 0), qt.basis(2, 1)
N4 = k0 * k0.dag()
N7 = k1 * k1.dag()
SX = k0 * k1.dag() + k1 * k0.dag()
SZ = N4 - N7
SM = k0 * k1.dag()
PSI0 = k0


def gen_ou(n, dt, var, tau, rng):
    """Exact OU discretization (analytic SDE solution)."""
    decay = np.exp(-dt / tau)
    noise = np.sqrt(var * (1 - decay**2)) * rng.standard_normal(n)
    x = np.empty(n)
    x[0] = np.sqrt(var) * rng.standard_normal()
    for i in range(1, n):
        x[i] = x[i - 1] * decay + noise[i]
    return x


def gen_noise(n, dt, sigma, amps, taus, rng):
    x = np.zeros(n)
    for A, tau in zip(amps, taus):
        if A > 0:
            x += gen_ou(n, dt, A * sigma**2, tau, rng)
    return x


def run_ensemble(sig4, sig7, amps4, taus4, amps7, taus7, n_real, H0):
    P4 = np.zeros((n_real, N_STEPS))
    P7 = np.zeros((n_real, N_STEPS))
    coh = np.zeros((n_real, N_STEPS))
    for i in range(n_real):
        rng = np.random.default_rng(42 + i)
        n4 = gen_noise(N_STEPS, DT, sig4, amps4, taus4, rng) * CM_TO_RADPS
        n7 = gen_noise(N_STEPS, DT, sig7, amps7, taus7, rng) * CM_TO_RADPS
        H_t = [H0, [N4, n4], [N7, n7]]
        r = qt.sesolve(H_t, PSI0, TLIST, e_ops=[N4, N7, SM])
        P4[i], P7[i], coh[i] = r.expect[0], r.expect[1], np.abs(r.expect[2])
    return dict(P4=P4.mean(0), P7=P7.mean(0), coh=coh.mean(0),
                P4_sem=P4.std(0)/np.sqrt(n_real),
                P7_sem=P7.std(0)/np.sqrt(n_real),
                coh_sem=coh.std(0)/np.sqrt(n_real))


def lindblad_2site():
    """Hakan-Strobl pure dephasing on the 2-site dimer, gamma=50 cm^-1."""
    J  = J_CM * CM_TO_RADPS
    dE = D_E_CM * CM_TO_RADPS
    g  = GAMMA_CM * CM_TO_RADPS
    H0 = J * SX + (dE / 2) * SZ
    c_ops = [np.sqrt(g) * N4, np.sqrt(g) * N7]
    r = qt.mesolve(H0, PSI0 * PSI0.dag(), TLIST, c_ops=c_ops,
                   e_ops=[N4, N7, SM])
    return dict(t=TLIST, P4=r.expect[0], P7=r.expect[1],
                coh=np.abs(r.expect[2]))


def lindblad_8site():
    """Full 8-site Craddock Hamiltonian, Hakan-Strobl gamma=50 cm^-1.
    Initial excitation on Trp4; observe all populations (matches notebook)."""
    H = qt.Qobj(H_CRAD * CM_TO_RADPS)
    g = GAMMA_CM * CM_TO_RADPS
    c_ops = [np.sqrt(g) * qt.basis(8, m) * qt.basis(8, m).dag()
             for m in range(8)]
    e_ops = [qt.basis(8, n) * qt.basis(8, n).dag() for n in range(8)]
    psi0 = qt.basis(8, I4)
    r = qt.mesolve(H, psi0, TLIST, c_ops=c_ops, e_ops=e_ops)
    return dict(t=TLIST, P=np.array(r.expect), donor=I4, target=I7)


def main():
    sig, tri = load_site_params()
    s4, s7 = sig[I4], sig[I7]
    a4, t4 = tri[I4]
    a7, t7 = tri[I7]
    J  = J_CM * CM_TO_RADPS
    dE = D_E_CM * CM_TO_RADPS
    H0 = J * SX + (dE / 2) * SZ

    print(f"Dimer Trp4(αW407)-Trp7(βW344)  J={J_CM} dE={D_E_CM} γ_L={GAMMA_CM} cm⁻¹")
    print(f"σ4={s4:.1f}  σ7={s7:.1f}  (slow traj)")
    print(f"Trp4  A={a4}  τ={t4}")
    print(f"Trp7  A={a7}  τ={t7}")
    print(f"Grid {N_STEPS} steps dt={DT*1e3:.0f}fs T={T_MAX}ps  N_real={N_REAL}\n")

    # noiseless reference
    r0 = qt.sesolve(H0, PSI0, TLIST, e_ops=[N4, N7, SM])
    print("noiseless done")

    # ablation configs. Each site uses its own amps; 'only' configs zero the
    # other components but keep that site's amplitude for the kept component.
    def cfg(keep):
        ka4 = a4.copy(); ka7 = a7.copy()
        for k in (0, 1, 2):
            if k not in keep:
                ka4[k] = 0.0; ka7[k] = 0.0
        return ka4, ka7
    configs = [
        ('tau1',      cfg({0})),
        ('tau2',      cfg({1})),
        ('tau1_tau2', cfg({0, 1})),
        ('full',      cfg({0, 1, 2})),
        ('tau1_tau3', cfg({0, 2})),
        ('tau3',      cfg({2})),
    ]
    out = dict(t=TLIST,
               P4_noiseless=r0.expect[0], P7_noiseless=r0.expect[1],
               coh_noiseless=np.abs(r0.expect[2]),
               J_cm=J_CM, dE_cm=D_E_CM, sigma4=s4, sigma7=s7)
    for name, (ka4, ka7) in configs:
        print(f"  MC {name:10s} ...", end=' ', flush=True)
        r = run_ensemble(s4, s7, ka4, t4, ka7, t7, N_REAL, H0)
        for k, v in r.items():
            out[f'{name}_{k}'] = v
        print(f"P7(2ps)={r['P7'][-1]:.3f}")
    np.savez_compressed(PAPER_DIR / 'exciton_mc.npz', **out)
    # conservation + noiseless checks
    for name, _ in configs:
        Psum = out[f'{name}_P4'] + out[f'{name}_P7']
        print(f"  {name}: max|P4+P7-1| = {np.max(np.abs(Psum-1)):.2e}")
    print(f"  noiseless P7 peak = {r0.expect[1].max():.4f} "
          f"(expect sin² peak {(1-0):.3f} at J-only)")

    print("\nLindblad 2-site (γ=50)...")
    lb2 = lindblad_2site()
    np.savez_compressed(PAPER_DIR / 'lindblad_2site.npz', **lb2)
    print(f"  Lindblad P7(2ps) = {lb2['P7'][-1]:.4f}")

    print("Lindblad 8-site Craddock (γ=50)...")
    lb8 = lindblad_8site()
    np.savez_compressed(PAPER_DIR / 'lindblad_8site.npz', **lb8)
    print(f"  Trp7 P(2ps) = {lb8['P'][I7][-1]:.4f}")

    print("\nsummary  P7(2ps):")
    print(f"  {'config':<12} {'MC':>8} {'Lindblad2':>10}")
    print(f"  {'full':<12} {out['full_P7'][-1]:>8.3f} {lb2['P7'][-1]:>10.3f}")
    print(f"  {'tau3':<12} {out['tau3_P7'][-1]:>8.3f}")
    print(f"  {'tau1':<12} {out['tau1_P7'][-1]:>8.3f}")


if __name__ == '__main__':
    main()
